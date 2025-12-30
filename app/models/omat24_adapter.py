"""
OMAT24 (FAIRChem/EquiformerV2) 模型适配器

GitHub: https://github.com/FAIR-Chem/fairchem
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class OMAT24Adapter(BaseModelAdapter):
    """
    OMAT24 (FAIRChem EquiformerV2) 模型适配器
    
    支持的模型:
    - omat24_prod_mp: eqV2_dens_86M_mp
    - omat24_prod: eqV2_86M_omat_mp_salex
    - omat24_prod_esen: esen_30m_oam
    - omat24_prod_esen_mp: esen_30m_mptrj
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "omat24_prod_mp": {
            "checkpoint": "eqV2_dens_86M_mp",
            "dispersion": True
        },
        "omat24_prod": {
            "checkpoint": "eqV2_86M_omat_mp_salex",
            "dispersion": True
        },
        "omat24_prod_esen": {
            "checkpoint": "esen_30m_oam",
            "dispersion": True
        },
        "omat24_prod_esen_mp": {
            "checkpoint": "esen_30m_mptrj",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 OMAT24 模型"""
        if self._is_loaded:
            return
        
        import torch
        from fairchem.core import OCPCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        checkpoint = config.get("checkpoint")
        
        self._calculator = OCPCalculator(
            model_name=checkpoint,
            local_cache="./models_cache",
            cpu=(self.device == "cpu")
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
        """获取 OMAT24 Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
