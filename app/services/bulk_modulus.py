"""
体积模量计算服务
"""
from typing import Optional, Callable, List, Tuple
from pathlib import Path
import uuid
import numpy as np
from scipy.optimize import curve_fit

from ase import Atoms
from ase.io import write

from ..models.base import BaseModelAdapter
from ..schemas.task import BulkModulusRequest
from ..schemas.result import BulkModulusResult
from ..core.ase_utils import get_volume_scaled_structures


class BulkModulusService:
    """体积模量计算服务"""
    
    def __init__(self, results_dir: Path):
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    @staticmethod
    def birch_murnaghan_eos(V, E0, V0, B0, B0_prime):
        """
        Birch-Murnaghan 状态方程
        
        Args:
            V: 体积
            E0: 平衡能量
            V0: 平衡体积
            B0: 体积模量
            B0_prime: 体积模量对压力的导数
            
        Returns:
            能量
        """
        eta = (V0 / V) ** (2/3)
        E = E0 + (9 * V0 * B0 / 16) * (
            (eta - 1) ** 3 * B0_prime + 
            (eta - 1) ** 2 * (6 - 4 * eta)
        )
        return E
    
    @staticmethod
    def murnaghan_eos(V, E0, V0, B0, B0_prime):
        """
        Murnaghan 状态方程
        """
        E = E0 + B0 * V / B0_prime * (
            ((V0 / V) ** B0_prime) / (B0_prime - 1) + 1
        ) - V0 * B0 / (B0_prime - 1)
        return E
    
    def run(
        self,
        atoms: Atoms,
        model: BaseModelAdapter,
        request: BulkModulusRequest,
        task_id: Optional[str] = None,
        progress_callback: Optional[Callable[[float], None]] = None
    ) -> BulkModulusResult:
        """
        计算体积模量
        
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
        
        # 设置计算器
        atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
        
        # 预优化（可选）
        if request.pre_optimize:
            from ase.optimize import BFGS
            from ase.filters import FrechetCellFilter
            opt = BFGS(FrechetCellFilter(atoms))
            opt.run(fmax=0.01, steps=200)
        
        # 生成不同体积的结构
        structures = get_volume_scaled_structures(
            atoms,
            volume_range=request.volume_range,
            num_points=request.num_points
        )
        
        # 计算每个体积对应的能量
        volumes = []
        energies = []
        
        for i, (vol, scaled_atoms) in enumerate(structures):
            scaled_atoms.calc = model.get_calculator(enable_d3=request.enable_d3)
            energy = scaled_atoms.get_potential_energy()
            volumes.append(vol)
            energies.append(energy)
            
            if progress_callback:
                progress = ((i + 1) / len(structures)) * 100
                progress_callback(progress)
        
        volumes = np.array(volumes)
        energies = np.array(energies)
        
        # 初始猜测
        V0_guess = atoms.get_volume()
        E0_guess = min(energies)
        B0_guess = 10.0  # GPa
        B0_prime_guess = 4.0
        
        # 拟合 Birch-Murnaghan 状态方程
        try:
            # 将体积模量转换为 eV/Å³
            # 1 GPa = 0.00624150913 eV/Å³
            GPa_to_eV_A3 = 0.00624150913
            
            popt, pcov = curve_fit(
                self.birch_murnaghan_eos,
                volumes,
                energies,
                p0=[E0_guess, V0_guess, B0_guess * GPa_to_eV_A3, B0_prime_guess],
                maxfev=10000
            )
            
            E0, V0, B0_eV, B0_prime = popt
            B0_GPa = B0_eV / GPa_to_eV_A3  # 转换回 GPa
            
            # 计算拟合质量
            E_fit = self.birch_murnaghan_eos(volumes, *popt)
            ss_res = np.sum((energies - E_fit) ** 2)
            ss_tot = np.sum((energies - np.mean(energies)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            residual = np.sqrt(ss_res / len(energies))
            
        except Exception as e:
            # 拟合失败时使用默认值
            E0 = E0_guess
            V0 = V0_guess
            B0_GPa = 0.0
            B0_prime = 4.0
            r_squared = 0.0
            residual = float('inf')
        
        # 构建 E-V 曲线数据
        ev_data = [
            {"volume": float(v), "energy": float(e)}
            for v, e in zip(volumes, energies)
        ]
        
        # 保存数据
        data_file = self.results_dir / f"{task_id}_ev_data.json"
        import json
        with open(data_file, 'w') as f:
            json.dump({
                "volumes": volumes.tolist(),
                "energies": energies.tolist(),
                "fit_params": {
                    "E0": float(E0),
                    "V0": float(V0),
                    "B0_GPa": float(B0_GPa),
                    "B0_prime": float(B0_prime)
                }
            }, f, indent=2)
        
        return BulkModulusResult(
            task_id=task_id,
            bulk_modulus=float(B0_GPa),
            bulk_modulus_derivative=float(B0_prime),
            equilibrium_volume=float(V0),
            equilibrium_energy=float(E0),
            r_squared=float(r_squared),
            residual=float(residual),
            ev_curve_data=ev_data,
            data_file=str(data_file)
        )
