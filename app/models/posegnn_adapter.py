"""
PosEGNN 模型适配器

PosEGNN: Position-based Equivariant Graph Neural Network
IBM Research 开发的基础模型

GitHub: https://github.com/IBM/materials/tree/main/models/pos_egnn
HuggingFace: https://huggingface.co/ibm-research/materials.pos-egnn
"""
from typing import Optional, Dict, Any
from .base import BaseModelAdapter

try:
    from ase.calculators.calculator import Calculator
except ImportError:
    Calculator = Any


class PosEGNNAdapter(BaseModelAdapter):
    """
    PosEGNN 模型适配器
    
    支持的模型:
    - posegnn_prod: pos-egnn.v1-6M
    """
    
    MODEL_CONFIG: Dict[str, Dict[str, Any]] = {
        "posegnn_prod": {
            "repo_id": "ibm-research/materials.pos-egnn",
            "filename": "pos-egnn.v1-6M.pt",
            "dispersion": True
        },
    }
    
    def load(self) -> None:
        """加载 PosEGNN 模型"""
        if self._is_loaded:
            return
        
        import torch
        from huggingface_hub import hf_hub_download
        
        # 延迟导入 posegnn
        try:
            from posegnn import PosEGNNCalculator
        except ImportError:
            # 尝试从 IBM materials 仓库导入
            import sys
            sys.path.append("./models_cache/pos_egnn")
            from load import PosEGNNCalculator
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        
        # 从 HuggingFace 下载模型权重
        model_path = hf_hub_download(
            repo_id=config.get("repo_id"),
            filename=config.get("filename")
        )
        
        self._calculator = PosEGNNCalculator(
            model_path=model_path,
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
        """获取 PosEGNN Calculator"""
        if not self._is_loaded:
            self.load()
        
        config = self.MODEL_CONFIG.get(self.model_key, {})
        needs_d3 = config.get("dispersion", True)
        
        if enable_d3 and needs_d3:
            return self._add_d3_correction(self._calculator)
        
        return self._calculator
