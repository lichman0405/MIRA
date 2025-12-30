"""
结构优化服务
"""
from typing import Optional, Callable, Dict, Any
from pathlib import Path
import uuid
import numpy as np

from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
# ASE 版本兼容性：ExpCellFilter 在 ASE >= 3.23 中移到了 ase.filters
try:
    from ase.filters import FrechetCellFilter, ExpCellFilter
except ImportError:
    from ase.constraints import ExpCellFilter
    from ase.filters import FrechetCellFilter
from ase.io import write
from ase.io.trajectory import Trajectory

from ..models.base import BaseModelAdapter
from ..schemas.task import OptimizationRequest
from ..schemas.result import OptimizationResult
from ..core.ase_utils import calculate_rmsd, get_structure_info


class OptimizationService:
    """结构优化服务"""
    
    OPTIMIZERS = {
        "BFGS": BFGS,
        "FIRE": FIRE,
        "LBFGS": LBFGS
    }
    
    FILTERS = {
        "frechet": FrechetCellFilter,
        "exp": ExpCellFilter
    }
    
    def __init__(self, results_dir: Path):
        """
        初始化优化服务
        
        Args:
            results_dir: 结果保存目录
        """
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def run(
        self,
        atoms: Atoms,
        model: BaseModelAdapter,
        request: OptimizationRequest,
        task_id: Optional[str] = None,
        progress_callback: Optional[Callable[[float], None]] = None
    ) -> OptimizationResult:
        """
        执行结构优化
        
        Args:
            atoms: ASE Atoms 对象
            model: 模型适配器
            request: 优化请求参数
            task_id: 任务 ID（可选）
            progress_callback: 进度回调函数
            
        Returns:
            优化结果
        """
        task_id = task_id or str(uuid.uuid4())
        
        # 复制初始结构
        initial_atoms = atoms.copy()
        initial_volume = atoms.get_volume()
        
        # 设置计算器
        atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
        
        # 计算初始能量
        initial_energy = atoms.get_potential_energy()
        
        # 应用晶胞过滤器
        if request.use_filter:
            atoms_to_optimize = FrechetCellFilter(atoms)
        else:
            atoms_to_optimize = atoms
        
        # 选择优化器
        OptimizerClass = self.OPTIMIZERS.get(request.optimizer, BFGS)
        
        # 轨迹文件
        traj_file = self.results_dir / f"{task_id}_opt.traj"
        log_file = self.results_dir / f"{task_id}_opt.log"
        
        optimizer = OptimizerClass(
            atoms_to_optimize,
            trajectory=str(traj_file),
            logfile=str(log_file)
        )
        
        # 进度追踪
        step_count = 0
        max_steps = request.max_steps
        
        def step_callback():
            nonlocal step_count
            step_count += 1
            if progress_callback:
                progress = min(100.0, (step_count / max_steps) * 100)
                progress_callback(progress)
        
        optimizer.attach(step_callback)
        
        # 运行优化
        try:
            converged = optimizer.run(fmax=request.fmax, steps=request.max_steps)
        except Exception as e:
            converged = False
        
        # 计算最终结果
        final_energy = atoms.get_potential_energy()
        final_volume = atoms.get_volume()
        final_forces = atoms.get_forces()
        final_fmax = np.max(np.linalg.norm(final_forces, axis=1))
        
        # 计算 RMSD
        rmsd = calculate_rmsd(initial_atoms, atoms)
        
        # 保存优化后的结构
        opt_file = self.results_dir / f"{task_id}_optimized.cif"
        write(str(opt_file), atoms)
        
        # 获取晶胞参数
        cell = atoms.get_cell()
        cell_params = {
            "a": float(cell.lengths()[0]),
            "b": float(cell.lengths()[1]),
            "c": float(cell.lengths()[2]),
            "alpha": float(cell.angles()[0]),
            "beta": float(cell.angles()[1]),
            "gamma": float(cell.angles()[2])
        }
        
        return OptimizationResult(
            task_id=task_id,
            initial_energy=float(initial_energy),
            final_energy=float(final_energy),
            energy_change=float(final_energy - initial_energy),
            energy_per_atom=float(final_energy / len(atoms)),
            rmsd=float(rmsd),
            initial_volume=float(initial_volume),
            final_volume=float(final_volume),
            volume_change_percent=float((final_volume - initial_volume) / initial_volume * 100),
            num_steps=step_count,
            converged=bool(converged if converged is not None else final_fmax < request.fmax),
            final_fmax=float(final_fmax),
            cell_parameters=cell_params,
            optimized_structure_file=str(opt_file),
            trajectory_file=str(traj_file)
        )
