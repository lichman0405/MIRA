"""
模型管理 API 路由
"""
from typing import List, Optional
from fastapi import APIRouter, HTTPException, status

from ...schemas.model import (
    ModelInfo, ModelFamily, ModelLoadRequest,
    AVAILABLE_MODELS
)
from ...dependencies import get_model_registry


router = APIRouter(prefix="/models", tags=["Models"])


@router.get("/", response_model=List[ModelInfo])
async def list_models(
    family: Optional[ModelFamily] = None
) -> List[ModelInfo]:
    """
    列出所有可用模型
    
    可选通过模型家族筛选
    """
    models = list(AVAILABLE_MODELS.values())
    
    if family:
        models = [m for m in models if m.family == family]
    
    return models


@router.get("/families", response_model=List[str])
async def list_model_families() -> List[str]:
    """列出所有模型家族"""
    return [f.value for f in ModelFamily]


@router.get("/{model_key}", response_model=ModelInfo)
async def get_model_info(model_key: str) -> ModelInfo:
    """获取指定模型的详细信息"""
    if model_key not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model '{model_key}' not found"
        )
    
    return AVAILABLE_MODELS[model_key]


@router.get("/{model_key}/status")
async def get_model_status(model_key: str):
    """获取模型加载状态"""
    if model_key not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model '{model_key}' not found"
        )
    
    registry = get_model_registry()
    is_loaded = registry.is_loaded(model_key)
    
    return {
        "model_key": model_key,
        "loaded": is_loaded
    }


@router.post("/{model_key}/load")
async def load_model(
    model_key: str,
    request: Optional[ModelLoadRequest] = None
):
    """
    预加载模型
    
    在提交任务前预先加载模型以减少首次运行延迟
    """
    if model_key not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model '{model_key}' not found"
        )
    
    registry = get_model_registry()
    
    try:
        device = request.device if request else "cuda"
        model = registry.get_or_load_model(model_key, device=device)
        
        return {
            "model_key": model_key,
            "status": "loaded",
            "device": device
        }
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to load model: {str(e)}"
        )


@router.post("/{model_key}/unload")
async def unload_model(model_key: str):
    """卸载模型以释放 GPU 内存"""
    if model_key not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Model '{model_key}' not found"
        )
    
    registry = get_model_registry()
    
    if not registry.is_loaded(model_key):
        return {
            "model_key": model_key,
            "status": "not_loaded"
        }
    
    registry.unload_model(model_key)
    
    return {
        "model_key": model_key,
        "status": "unloaded"
    }


@router.get("/loaded/list")
async def list_loaded_models():
    """列出当前已加载的模型"""
    registry = get_model_registry()
    loaded = registry.list_loaded_models()
    
    return {
        "loaded_models": loaded,
        "count": len(loaded)
    }
