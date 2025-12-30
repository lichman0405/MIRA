"""
MD 稳定性模拟服务
"""
from typing import Optional, Callable, Dict, Any, List
from pathlib import Path
import uuid
import numpy as np

from ase import Atoms, units
from ase.md.langevin import Langevin
from ase.md.nptberendsen import NPTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io import write
from ase.io.trajectory import Trajectory

from ..models.base import BaseModelAdapter
from ..schemas.task import StabilityRequest
from ..schemas.result import StabilityResult
from ..core.ase_utils import calculate_rmsd, get_coordination_numbers
from .optimization import OptimizationService


class StabilityService:
    """MD 稳定性模拟服务"""
    
    def __init__(self, results_dir: Path):
        """
        初始化服务
        
        Args:
            results_dir: 结果保存目录
        """
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def run(
        self,
        atoms: Atoms,
        model: BaseModelAdapter,
        request: StabilityRequest,
        task_id: Optional[str] = None,
        progress_callback: Optional[Callable[[float], None]] = None
    ) -> StabilityResult:
        """
        执行 MD 稳定性模拟
        
        Args:
            atoms: ASE Atoms 对象
            model: 模型适配器
            request: 模拟请求参数
            task_id: 任务 ID
            progress_callback: 进度回调
            
        Returns:
            模拟结果
        """
        task_id = task_id or str(uuid.uuid4())
        
        # 复制初始结构
        initial_atoms = atoms.copy()
        
        # 设置计算器
        atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
        
        # Step 1: 预优化（可选）
        if request.pre_optimize:
            from ..schemas.task import OptimizationRequest
            opt_service = OptimizationService(self.results_dir)
            opt_request = OptimizationRequest(
                structure_id="",
                model_key=request.model_key,
                fmax=request.fmax,
                max_steps=200,
                optimizer="BFGS",
                use_filter=True,
                enable_d3=request.enable_d3
            )
            # 简单优化
            from ase.optimize import BFGS
            from ase.filters import FrechetCellFilter
            opt = BFGS(FrechetCellFilter(atoms))
            opt.run(fmax=request.fmax, steps=200)
        
        # Step 2: 初始化速度
        MaxwellBoltzmannDistribution(atoms, temperature_K=request.nvt_temperature)
        
        # 数据记录
        temperatures = []
        volumes = []
        energies = []
        times = []
        
        total_steps = request.nvt_steps + request.npt_steps
        current_step = 0
        
        def record_data():
            nonlocal current_step
            current_step += 1
            temperatures.append(atoms.get_temperature())
            volumes.append(atoms.get_volume())
            energies.append(atoms.get_potential_energy())
            times.append(current_step * request.timestep)
            
            if progress_callback:
                progress = (current_step / total_steps) * 100
                progress_callback(min(progress, 100.0))
        
        # Step 3: NVT 平衡
        nvt_traj = self.results_dir / f"{task_id}_nvt.traj"
        nvt = Langevin(
            atoms,
            timestep=request.timestep * units.fs,
            temperature_K=request.nvt_temperature,
            friction=0.01 / units.fs,
            trajectory=str(nvt_traj)
        )
        nvt.attach(record_data, interval=request.trajectory_interval)
        nvt.run(request.nvt_steps)
        
        # Step 4: NPT 生产运行
        npt_traj = self.results_dir / f"{task_id}_npt.traj"
        
        # 使用 NPTBerendsen（ASE 3.27+ 支持）
        npt = NPTBerendsen(
            atoms,
            timestep=request.timestep * units.fs,
            temperature_K=request.npt_temperature,
            pressure_au=request.npt_pressure * units.bar,
            taut=100 * units.fs,
            taup=1000 * units.fs,
            compressibility_au=4.57e-5 / units.bar,
            trajectory=str(npt_traj)
        )
        npt.attach(record_data, interval=request.trajectory_interval)
        npt.run(request.npt_steps)
        
        # 计算分析结果
        final_rmsd = calculate_rmsd(initial_atoms, atoms)
        
        # 计算 RMSD 随时间的演化
        traj = Trajectory(str(npt_traj))
        rmsd_evolution = []
        for frame in traj:
            rmsd_evolution.append(calculate_rmsd(initial_atoms, frame))
        max_rmsd = max(rmsd_evolution) if rmsd_evolution else final_rmsd
        
        # 配位数分析
        initial_coord = get_coordination_numbers(initial_atoms)
        final_coord = get_coordination_numbers(atoms)
        
        coord_analysis = {
            "initial": initial_coord,
            "final": final_coord
        }
        
        # 保存最终结构
        final_file = self.results_dir / f"{task_id}_final.cif"
        write(str(final_file), atoms)
        
        # 体积漂移
        initial_volume = initial_atoms.get_volume()
        final_volume = atoms.get_volume()
        volume_drift = (final_volume - initial_volume) / initial_volume * 100
        
        return StabilityResult(
            task_id=task_id,
            average_temperature=float(np.mean(temperatures)),
            temperature_std=float(np.std(temperatures)),
            average_pressure=float(request.npt_pressure),
            pressure_std=0.0,  # 需要从 MD 中提取
            final_rmsd=float(final_rmsd),
            max_rmsd=float(max_rmsd),
            volume_drift_percent=float(volume_drift),
            coordination_analysis=coord_analysis,
            average_energy=float(np.mean(energies)),
            energy_std=float(np.std(energies)),
            total_steps=total_steps,
            simulation_time_ps=float(total_steps * request.timestep / 1000),
            trajectory_file=str(npt_traj),
            final_structure_file=str(final_file)
        )
