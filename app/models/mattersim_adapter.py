"""
MatterSim 模型适配器

MatterSim: Microsoft Materials Simulation
GitHub: https://github.com/microsoft/mattersim
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class MatterSimAdapter(BaseModelAdapter):
    """
    MatterSim 模型适配器
    
    支持的模型:
    - mattersim_prod: MatterSim-v1.0.0-5M
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "mattersim_prod": {
            "model_name": "MatterSim-v1.0.0-5M",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 MatterSim 模型"""
        if self._is_loaded:
            return
        
        import torch
        from mattersim.forcefield import MatterSimCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        
        self._calculator = MatterSimCalculator(
            device=self.device
        )
        
        self._is_loaded = True
        self._memory_usage_mb = self._estimate_memory_usage()
    
    def unload(self) -> None:
        """卸载模型"""
        if not self._is_loaded:
            return
        
        import torch
        
        del self._calculator
        self._calculator = None
        self._is_loaded = False
        self._memory_usage_mb = None
        
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    
    def get_calculator(self, enable_d3: bool = True) -> Calculator:
        """获取 MatterSim Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
