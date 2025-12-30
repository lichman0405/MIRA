"""
GRACE 模型适配器

GRACE: Graph-based Atomic Cluster Expansion
相关: https://github.com/ACEsuit/mace
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class GRACEAdapter(BaseModelAdapter):
    """
    GRACE 模型适配器
    
    支持的模型:
    - grace_prod: GRACE-2L-MP-r6
    - grace_prod_oam: GRACE-2L-OAM
    - grace_prod_omat: GRACE-2L-OMAT
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "grace_prod": {
            "model_name": "GRACE-2L-MP-r6",
            "dispersion": True
        },
        "grace_prod_oam": {
            "model_name": "GRACE-2L-OAM",
            "dispersion": True
        },
        "grace_prod_omat": {
            "model_name": "GRACE-2L-OMAT",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 GRACE 模型"""
        if self._is_loaded:
            return
        
        import torch
        from tensorpotential.calculator import TPCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        model_name = config.get("model_name")
        
        self._calculator = TPCalculator(
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
        """获取 GRACE Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
