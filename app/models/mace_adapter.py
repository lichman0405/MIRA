"""
MACE 模型适配器

MACE: Multi-ACE (Atomic Cluster Expansion)
GitHub: https://github.com/ACEsuit/mace
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class MACEAdapter(BaseModelAdapter):
    """
    MACE 模型适配器
    
    支持的模型:
    - mace_prod_b3: MACE-MP-0b3-medium (官方预训练)
    - mace_prod: MACE-MPA-0-medium
    - mace_prod_omat: MACE-OMAT-0-medium
    - mace_prod_matpes: MACE-MatPES-r2scan
    - mace_prod_mof: MACE-MOF-v1
    """
    
    # 模型配置映射
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "mace_prod_b3": {
            "method": "mace_mp",
            "model": "medium",
            "dispersion": False  # MACE-MP 已包含色散
        },
        "mace_prod": {
            "method": "mace_mp",
            "model": "medium",
            "dispersion": False
        },
        "mace_prod_omat": {
            "method": "huggingface",
            "repo_id": "mace-ml/MACE-OMAT-0-medium",
            "dispersion": True
        },
        "mace_prod_matpes": {
            "method": "huggingface",
            "repo_id": "mace-ml/MACE-MatPES-r2scan",
            "dispersion": True
        },
        "mace_prod_mof": {
            "method": "local",
            "model_path": None,  # 需要用户提供
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 MACE 模型"""
        if self._is_loaded:
            return
            
        import torch
        from mace.calculators import mace_mp, MACECalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        method = config.get("method", "mace_mp")
        
        # 确定 dtype
        default_dtype = "float32" if self.dtype == "float32" else "float64"
        
        if method == "mace_mp":
            # 使用官方 mace_mp 函数自动下载
            self._calculator = mace_mp(
                model=config.get("model", "medium"),
                device=self.device,
                default_dtype=default_dtype
            )
        elif method == "huggingface":
            # 从 HuggingFace 下载
            from huggingface_hub import hf_hub_download
            repo_id = config.get("repo_id")
            model_path = hf_hub_download(
                repo_id=repo_id,
                filename="model.pt"
            )
            self._calculator = MACECalculator(
                model_paths=model_path,
                device=self.device,
                default_dtype=default_dtype
            )
        elif method == "local":
            # 本地模型文件
            model_path = config.get("model_path")
            if model_path is None:
                raise ValueError(f"Model path not configured for {self.model_key}")
            self._calculator = MACECalculator(
                model_paths=model_path,
                device=self.device,
                default_dtype=default_dtype
            )
        
        self._is_loaded = True
        self._memory_usage_mb = self._estimate_memory_usage()
    
    def unload(self) -> None:
        """卸载模型释放内存"""
        if not self._is_loaded:
            return
            
        import torch
        
        del self._calculator
        self._calculator = None
        self._is_loaded = False
        self._memory_usage_mb = None
        
        # 清理 GPU 缓存
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    
    def get_calculator(self, enable_d3: bool = True) -> Calculator:
        """
        获取 MACE Calculator
        
        Args:
            enable_d3: 是否启用 D3 校正（注意：部分 MACE 模型已包含色散）
        """
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
