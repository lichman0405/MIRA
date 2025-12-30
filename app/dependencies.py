"""
依赖注入
"""
from functools import lru_cache
from pathlib import Path

from .config import Settings, get_settings
from .models.registry import ModelRegistry
from .services.structure_service import StructureService
from .services.task_service import TaskService


# 全局单例
_model_registry: ModelRegistry = None
_structure_service: StructureService = None
_task_service: TaskService = None


def get_model_registry() -> ModelRegistry:
    """获取模型注册表单例"""
    global _model_registry
    if _model_registry is None:
        settings = get_settings()
        _model_registry = ModelRegistry(default_device=settings.DEFAULT_DEVICE)
    return _model_registry


def get_structure_service() -> StructureService:
    """获取结构服务单例"""
    global _structure_service
    if _structure_service is None:
        settings = get_settings()
        _structure_service = StructureService(
            structures_dir=Path(settings.STRUCTURES_DIR)
        )
    return _structure_service


def get_task_service() -> TaskService:
    """获取任务服务单例"""
    global _task_service
    if _task_service is None:
        settings = get_settings()
        _task_service = TaskService(
            model_registry=get_model_registry(),
            structures_dir=Path(settings.STRUCTURES_DIR),
            results_dir=Path(settings.RESULTS_DIR),
            max_workers=settings.MAX_CONCURRENT_TASKS
        )
    return _task_service


async def startup_handler():
    """应用启动时的初始化"""
    settings = get_settings()
    
    # 确保目录存在
    Path(settings.STRUCTURES_DIR).mkdir(parents=True, exist_ok=True)
    Path(settings.RESULTS_DIR).mkdir(parents=True, exist_ok=True)
    
    # 初始化服务
    get_model_registry()
    get_structure_service()
    get_task_service()


async def shutdown_handler():
    """应用关闭时的清理"""
    global _model_registry
    
    if _model_registry is not None:
        # 卸载所有模型
        for model_key in _model_registry.list_loaded_models():
            _model_registry.unload_model(model_key)
