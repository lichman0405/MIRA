"""
SevenNet (7net) 模型适配器

SevenNet: 7-layer Neural Network Potential
GitHub: https://github.com/MDIL-SNU/SevenNet
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class SevenNetAdapter(BaseModelAdapter):
    """
    SevenNet 模型适配器
    
    支持的模型:
    - sevennet_prod: 7net-0
    - sevennet_prod_l3i5: 7net-l3i5
    - sevennet_prod_ompa: 7net-mf-ompa
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "sevennet_prod": {
            "model_name": "7net-0",
            "dispersion": True
        },
        "sevennet_prod_l3i5": {
            "model_name": "7net-l3i5",
            "dispersion": True
        },
        "sevennet_prod_ompa": {
            "model_name": "7net-mf-ompa",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 SevenNet 模型"""
        if self._is_loaded:
            return
        
        import torch
        from sevenn.sevennet_calculator import SevenNetCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        model_name = config.get("model_name")
        
        self._calculator = SevenNetCalculator(
            model=model_name,
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
        """获取 SevenNet Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
