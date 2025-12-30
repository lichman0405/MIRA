"""
热容计算服务（基于声子计算）
"""
from typing import Optional, Callable, List
from pathlib import Path
import uuid
import numpy as np

from ase import Atoms
from ase.io import write

from ..models.base import BaseModelAdapter
from ..schemas.task import HeatCapacityRequest
from ..schemas.result import HeatCapacityResult
from ..core.ase_utils import create_supercell


class HeatCapacityService:
    """热容计算服务"""
    
    def __init__(self, results_dir: Path):
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def run(
        self,
        atoms: Atoms,
        model: BaseModelAdapter,
        request: HeatCapacityRequest,
        task_id: Optional[str] = None,
        progress_callback: Optional[Callable[[float], None]] = None
    ) -> HeatCapacityResult:
        """
        计算热容
        
        使用 Phonopy 进行声子计算，然后计算热力学性质
        
        Args:
            atoms: ASE Atoms 对象
            model: 模型适配器
            request: 计算请求参数
            task_id: 任务 ID
            progress_callback: 进度回调
            
        Returns:
            计算结果
        """
        task_id = task_id or str(uuid.uuid4())
        
        try:
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
        except ImportError:
            raise ImportError("Phonopy is required: pip install phonopy")
        
        # 设置计算器
        atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
        
        # 预优化
        if request.pre_optimize:
            from ase.optimize import BFGS
            from ase.filters import FrechetCellFilter
            opt = BFGS(FrechetCellFilter(atoms))
            opt.run(fmax=0.01, steps=200)
        
        # 转换为 Phonopy 格式
        phonopy_atoms = PhonopyAtoms(
            symbols=atoms.get_chemical_symbols(),
            cell=atoms.get_cell(),
            scaled_positions=atoms.get_scaled_positions()
        )
        
        # 创建 Phonopy 对象
        phonon = Phonopy(
            phonopy_atoms,
            supercell_matrix=np.diag(request.supercell),
            primitive_matrix='auto'
        )
        
        # 生成位移
        phonon.generate_displacements(distance=request.displacement)
        supercells = phonon.supercells_with_displacements
        
        if progress_callback:
            progress_callback(10.0)
        
        # 计算力
        forces = []
        for i, scell in enumerate(supercells):
            # 转换回 ASE Atoms
            sc_atoms = Atoms(
                symbols=scell.symbols,
                cell=scell.cell,
                scaled_positions=scell.scaled_positions,
                pbc=True
            )
            sc_atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
            f = sc_atoms.get_forces()
            forces.append(f)
            
            if progress_callback:
                progress = 10.0 + (i + 1) / len(supercells) * 60.0
                progress_callback(progress)
        
        # 设置力
        phonon.forces = forces
        
        # 产生力常数
        phonon.produce_force_constants()
        
        if progress_callback:
            progress_callback(75.0)
        
        # 设置 mesh
        phonon.run_mesh(request.mesh)
        
        # 计算热力学性质
        t_min, t_max = request.temperature_range
        temperatures = np.linspace(t_min, t_max, request.temperature_points)
        
        phonon.run_thermal_properties(
            t_min=t_min,
            t_max=t_max,
            t_step=(t_max - t_min) / (request.temperature_points - 1)
        )
        
        thermal = phonon.get_thermal_properties_dict()
        
        if progress_callback:
            progress_callback(90.0)
        
        # 提取热容数据
        cv_data = thermal.get('heat_capacity', [])
        entropy_data = thermal.get('entropy', [])
        free_energy_data = thermal.get('free_energy', [])
        temp_data = thermal.get('temperatures', temperatures)
        
        # 找到 300K 的值
        idx_300k = np.argmin(np.abs(np.array(temp_data) - 300))
        cv_300k = cv_data[idx_300k] if len(cv_data) > idx_300k else 0.0
        entropy_300k = entropy_data[idx_300k] if len(entropy_data) > idx_300k else 0.0
        free_energy_300k = free_energy_data[idx_300k] if len(free_energy_data) > idx_300k else 0.0
        
        # 获取声子频率
        frequencies = phonon.get_mesh_dict()['frequencies'].flatten()
        
        # 检查虚频
        imaginary_mask = frequencies < 0
        has_imaginary = np.any(imaginary_mask)
        num_imaginary = int(np.sum(imaginary_mask))
        
        # 构建 T-Cv 数据
        temperature_cv_data = [
            {"temperature": float(t), "cv": float(cv)}
            for t, cv in zip(temp_data, cv_data)
        ]
        
        # 保存数据
        data_file = self.results_dir / f"{task_id}_thermal.json"
        import json
        with open(data_file, 'w') as f:
            json.dump({
                "temperatures": list(temp_data),
                "heat_capacity": list(cv_data),
                "entropy": list(entropy_data),
                "free_energy": list(free_energy_data),
                "frequencies_THz": frequencies.tolist()
            }, f, indent=2)
        
        if progress_callback:
            progress_callback(100.0)
        
        return HeatCapacityResult(
            task_id=task_id,
            cv_300k=float(cv_300k),
            entropy_300k=float(entropy_300k),
            free_energy_300k=float(free_energy_300k),
            internal_energy_300k=float(free_energy_300k + 300 * entropy_300k / 1000),
            temperature_cv_data=temperature_cv_data,
            phonon_frequencies=frequencies[:100].tolist(),  # 前100个频率
            has_imaginary_modes=bool(has_imaginary),
            num_imaginary_modes=num_imaginary,
            data_file=str(data_file)
        )
