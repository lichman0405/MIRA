"""
ML 力场模型适配器基类
所有模型适配器都继承此基类，提供统一的接口
"""
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any
import numpy as np

# 使用延迟导入避免在未安装 ASE 时报错
try:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
except ImportError:
    Atoms = Any
    Calculator = Any


class BaseModelAdapter(ABC):
    """
    所有 ML 力场模型的基类
    
    子类需要实现:
    - load(): 加载模型到内存
    - unload(): 卸载模型释放内存
    - get_calculator(): 返回 ASE Calculator
    """
    
    def __init__(
        self, 
        model_key: str, 
        device: str = "cuda", 
        dtype: str = "float32"
    ):
        """
        初始化适配器
        
        Args:
            model_key: 模型唯一标识符
            device: 运行设备 (cuda/cpu)
            dtype: 数据类型 (float32/float64)
        """
        self.model_key = model_key
        self.device = device
        self.dtype = dtype
        self._calculator: Optional[Calculator] = None
        self._is_loaded = False
        self._memory_usage_mb: Optional[float] = None
    
    @abstractmethod
    def load(self) -> None:
        """加载模型到内存"""
        pass
    
    @abstractmethod
    def unload(self) -> None:
        """卸载模型释放内存"""
        pass
    
    @abstractmethod
    def get_calculator(self, enable_d3: bool = True) -> Calculator:
        """
        返回 ASE Calculator
        
        Args:
            enable_d3: 是否启用 D3 色散校正
            
        Returns:
            ASE Calculator 对象
        """
        pass
    
    @property
    def is_loaded(self) -> bool:
        """模型是否已加载"""
        return self._is_loaded
    
    @property
    def memory_usage_mb(self) -> Optional[float]:
        """内存使用量 (MB)"""
        return self._memory_usage_mb
    
    def calculate_energy(self, atoms: Atoms, enable_d3: bool = True) -> float:
        """
        计算结构能量
        
        Args:
            atoms: ASE Atoms 对象
            enable_d3: 是否启用 D3 校正
            
        Returns:
            能量值 (eV)
        """
        atoms.calc = self.get_calculator(enable_d3=enable_d3)
        return atoms.get_potential_energy()
    
    def calculate_forces(self, atoms: Atoms, enable_d3: bool = True) -> np.ndarray:
        """
        计算原子受力
        
        Args:
            atoms: ASE Atoms 对象
            enable_d3: 是否启用 D3 校正
            
        Returns:
            力数组 (N, 3) (eV/Å)
        """
        atoms.calc = self.get_calculator(enable_d3=enable_d3)
        return atoms.get_forces()
    
    def calculate_stress(self, atoms: Atoms, enable_d3: bool = True) -> np.ndarray:
        """
        计算应力张量
        
        Args:
            atoms: ASE Atoms 对象
            enable_d3: 是否启用 D3 校正
            
        Returns:
            应力张量 (6,) Voigt 表示 (eV/Å³)
        """
        atoms.calc = self.get_calculator(enable_d3=enable_d3)
        return atoms.get_stress()
    
    def single_point(
        self, 
        atoms: Atoms, 
        enable_d3: bool = True
    ) -> Dict[str, Any]:
        """
        单点计算：同时获取能量、力和应力
        
        Args:
            atoms: ASE Atoms 对象
            enable_d3: 是否启用 D3 校正
            
        Returns:
            包含 energy, forces, stress 的字典
        """
        atoms.calc = self.get_calculator(enable_d3=enable_d3)
        return {
            "energy": atoms.get_potential_energy(),
            "forces": atoms.get_forces(),
            "stress": atoms.get_stress()
        }
    
    def _add_d3_correction(self, calculator: Calculator) -> Calculator:
        """
        添加 D3 色散校正
        
        Args:
            calculator: 原始计算器
            
        Returns:
            带 D3 校正的组合计算器
        """
        try:
            from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
            from ase.calculators.mixing import SumCalculator
            
            d3_calc = TorchDFTD3Calculator(
                device=self.device,
                damping="bj"
            )
            return SumCalculator([calculator, d3_calc])
        except ImportError:
            # 如果 torch-dftd 未安装，返回原始计算器
            import warnings
            warnings.warn("torch-dftd not installed, D3 correction disabled")
            return calculator
    
    def _estimate_memory_usage(self) -> float:
        """估算 GPU 内存使用量 (MB)"""
        try:
            import torch
            if torch.cuda.is_available() and self.device == "cuda":
                return torch.cuda.memory_allocated() / 1024 / 1024
        except ImportError:
            pass
        return 0.0
    
    def __repr__(self) -> str:
        status = "loaded" if self._is_loaded else "not loaded"
        return f"<{self.__class__.__name__}(key={self.model_key}, device={self.device}, {status})>"
