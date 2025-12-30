"""
模型注册表 - 管理所有模型适配器的加载和访问
"""
from typing import Dict, List, Optional, Type
from .base import BaseModelAdapter
from .mace_adapter import MACEAdapter
from .orb_adapter import ORBAdapter
from .omat24_adapter import OMAT24Adapter
from .grace_adapter import GRACEAdapter
from .mattersim_adapter import MatterSimAdapter
from .sevennet_adapter import SevenNetAdapter
from .posegnn_adapter import PosEGNNAdapter
from .matgl_adapter import MatGLAdapter
from ..schemas.model import ModelFamily, ModelInfo, AVAILABLE_MODELS


class ModelRegistry:
    """
    模型注册表
    
    管理所有 ML 力场模型的加载、卸载和访问
    """
    
    # 模型系列到适配器类的映射
    ADAPTER_CLASSES: Dict[ModelFamily, Type[BaseModelAdapter]] = {
        ModelFamily.MACE: MACEAdapter,
        ModelFamily.ORB: ORBAdapter,
        ModelFamily.OMAT24: OMAT24Adapter,
        ModelFamily.GRACE: GRACEAdapter,
        ModelFamily.MATTERSIM: MatterSimAdapter,
        ModelFamily.SEVENNET: SevenNetAdapter,
        ModelFamily.POSEGNN: PosEGNNAdapter,
        ModelFamily.MATGL: MatGLAdapter,
    }

    
    def __init__(self):
        """初始化模型注册表"""
        self._loaded_models: Dict[str, BaseModelAdapter] = {}
        self._default_device = "cuda"
        self._default_dtype = "float32"
    
    def set_defaults(self, device: str = "cuda", dtype: str = "float32") -> None:
        """设置默认设备和数据类型"""
        self._default_device = device
        self._default_dtype = dtype
    
    def list_available_models(self) -> List[ModelInfo]:
        """列出所有可用模型"""
        models = []
        for key, config in AVAILABLE_MODELS.items():
            is_loaded = key in self._loaded_models
            adapter = self._loaded_models.get(key)
            models.append(ModelInfo(
                key=key,
                name=config["name"],
                family=config["family"],
                description=config.get("description"),
                supports_d3=config.get("supports_d3", True),
                is_loaded=is_loaded,
                memory_usage_mb=adapter.memory_usage_mb if adapter else None,
                device=adapter.device if adapter else None
            ))
        return models
    
    def get_model_info(self, model_key: str) -> Optional[ModelInfo]:
        """获取指定模型的信息"""
        if model_key not in AVAILABLE_MODELS:
            return None
        
        config = AVAILABLE_MODELS[model_key]
        is_loaded = model_key in self._loaded_models
        adapter = self._loaded_models.get(model_key)
        
        return ModelInfo(
            key=model_key,
            name=config["name"],
            family=config["family"],
            description=config.get("description"),
            supports_d3=config.get("supports_d3", True),
            is_loaded=is_loaded,
            memory_usage_mb=adapter.memory_usage_mb if adapter else None,
            device=adapter.device if adapter else None
        )
    
    def load_model(
        self, 
        model_key: str, 
        device: Optional[str] = None,
        dtype: Optional[str] = None
    ) -> BaseModelAdapter:
        """
        加载模型
        
        Args:
            model_key: 模型标识符
            device: 设备 (cuda/cpu)，默认使用注册表默认值
            dtype: 数据类型，默认使用注册表默认值
            
        Returns:
            加载后的模型适配器
        """
        # 如果已加载，直接返回
        if model_key in self._loaded_models:
            return self._loaded_models[model_key]
        
        # 检查模型是否存在
        if model_key not in AVAILABLE_MODELS:
            raise ValueError(f"Unknown model: {model_key}")
        
        config = AVAILABLE_MODELS[model_key]
        family = config["family"]
        
        # 获取适配器类
        adapter_class = self.ADAPTER_CLASSES.get(family)
        if adapter_class is None:
            raise ValueError(f"No adapter for model family: {family}")
        
        # 创建适配器
        adapter = adapter_class(
            model_key=model_key,
            device=device or self._default_device,
            dtype=dtype or self._default_dtype
        )
        
        # 加载模型
        adapter.load()
        
        # 注册到已加载模型
        self._loaded_models[model_key] = adapter
        
        return adapter
    
    def unload_model(self, model_key: str) -> bool:
        """
        卸载模型
        
        Args:
            model_key: 模型标识符
            
        Returns:
            是否成功卸载
        """
        if model_key not in self._loaded_models:
            return False
        
        adapter = self._loaded_models[model_key]
        adapter.unload()
        del self._loaded_models[model_key]
        
        return True
    
    def unload_all(self) -> None:
        """卸载所有已加载的模型"""
        for model_key in list(self._loaded_models.keys()):
            self.unload_model(model_key)
    
    def get_model(self, model_key: str) -> Optional[BaseModelAdapter]:
        """
        获取已加载的模型
        
        Args:
            model_key: 模型标识符
            
        Returns:
            模型适配器，如果未加载则返回 None
        """
        return self._loaded_models.get(model_key)
    
    def get_or_load_model(
        self, 
        model_key: str,
        device: Optional[str] = None,
        dtype: Optional[str] = None
    ) -> BaseModelAdapter:
        """
        获取模型，如果未加载则先加载
        
        Args:
            model_key: 模型标识符
            device: 设备
            dtype: 数据类型
            
        Returns:
            模型适配器
        """
        if model_key in self._loaded_models:
            return self._loaded_models[model_key]
        return self.load_model(model_key, device, dtype)
    
    @property
    def loaded_models(self) -> List[str]:
        """已加载的模型列表"""
        return list(self._loaded_models.keys())
    
    @property
    def total_memory_mb(self) -> float:
        """所有已加载模型的总内存使用量"""
        total = 0.0
        for adapter in self._loaded_models.values():
            if adapter.memory_usage_mb:
                total += adapter.memory_usage_mb
        return total
    
    def __contains__(self, model_key: str) -> bool:
        """检查模型是否已加载"""
        return model_key in self._loaded_models
    
    def __len__(self) -> int:
        """已加载模型数量"""
        return len(self._loaded_models)


# 全局模型注册表实例
model_registry = ModelRegistry()
