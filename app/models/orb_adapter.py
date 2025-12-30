"""
ORB 模型适配器

ORB: Orbital Materials Foundation Model
GitHub: https://github.com/orbital-materials/orb-models
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class ORBAdapter(BaseModelAdapter):
    """
    ORB 模型适配器
    
    支持的模型:
    - orb_prod_mp: orb-mptraj-only-v2
    - orb_prod: orb-d3-v2
    - orb3: orb-v3-direct-20-omat
    - orb_prod_v3: orb-v3-conservative-inf-omat
    - orb_prod_v3_mp: orb-v3-conservative-inf-mpa
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "orb_prod_mp": {
            "model_name": "orb-mptraj-only-v2",
            "dispersion": True
        },
        "orb_prod": {
            "model_name": "orb-d3-v2",
            "dispersion": False  # 已包含 D3
        },
        "orb3": {
            "model_name": "orb-v3-direct-20-omat",
            "dispersion": True
        },
        "orb_prod_v3": {
            "model_name": "orb-v3-conservative-inf-omat",
            "dispersion": True
        },
        "orb_prod_v3_mp": {
            "model_name": "orb-v3-conservative-inf-mpa",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 ORB 模型"""
        if self._is_loaded:
            return
        
        import torch
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        model_name = config.get("model_name")
        
        # 加载预训练模型
        orbff = pretrained.orb_v2(weights_path=model_name)
        
        self._calculator = ORBCalculator(
            orbff,
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
        """获取 ORB Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
