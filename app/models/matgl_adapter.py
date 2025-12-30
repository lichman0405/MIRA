"""
MatGL 模型适配器 (M3GNet, CHGNet)
"""
from typing import Optional, Any
from ase.calculators.calculator import Calculator

from .base import BaseModelAdapter
from ..schemas.model import AVAILABLE_MODELS


class MatGLAdapter(BaseModelAdapter):
    """MatGL 模型适配器"""
    
    def __init__(
        self,
        model_key: str,
        device: str = "cuda",
        **kwargs
    ):
        super().__init__(model_key, device, **kwargs)
        self._calculator = None
        self._potential = None
    
    def load(self) -> None:
        """加载 MatGL 模型"""
        if self._loaded:
            return
        
        try:
            import matgl
            from matgl.ext.ase import M3GNetCalculator, CHGNetCalculator
        except ImportError:
            raise ImportError("MatGL is required: pip install matgl")
        
        model_info = AVAILABLE_MODELS.get(self.model_key, {})
        model_path = model_info.model_path if hasattr(model_info, 'model_path') else None
        
        # 根据模型类型加载
        if "m3gnet" in self.model_key.lower():
            # 加载 M3GNet
            if model_path:
                pot = matgl.load_model(model_path)
            else:
                pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            self._potential = pot
            self._calculator = M3GNetCalculator(potential=pot)
        
        elif "chgnet" in self.model_key.lower():
            # 加载 CHGNet
            try:
                from chgnet.model import CHGNet
                from chgnet.model.dynamics import CHGNetCalculator as CHGNetCalc
                
                model = CHGNet.load()
                self._potential = model
                self._calculator = CHGNetCalc(model=model, use_device=self.device)
            except ImportError:
                # 回退到 MatGL 的 CHGNet
                if model_path:
                    pot = matgl.load_model(model_path)
                else:
                    pot = matgl.load_model("CHGNet-v0.3.0")
                self._potential = pot
                self._calculator = CHGNetCalculator(potential=pot)
        
        else:
            # 默认使用 M3GNet
            pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            self._potential = pot
            self._calculator = M3GNetCalculator(potential=pot)
        
        self._loaded = True
    
    def unload(self) -> None:
        """卸载模型"""
        if not self._loaded:
            return
        
        import gc
        import torch
        
        self._calculator = None
        self._potential = None
        
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        self._loaded = False
    
    def get_calculator(self, enable_d3: bool = True) -> Calculator:
        """获取 ASE Calculator"""
        if not self._loaded:
            self.load()
        
        if enable_d3:
            return self._get_d3_corrected_calculator(self._calculator)
        
        return self._calculator
