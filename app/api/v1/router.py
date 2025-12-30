"""
API v1 路由聚合
"""
from fastapi import APIRouter

from .models import router as models_router
from .structures import router as structures_router
from .tasks import router as tasks_router
from .results import router as results_router


# 创建 v1 路由
api_router = APIRouter(prefix="/v1")

# 包含所有子路由
api_router.include_router(models_router)
api_router.include_router(structures_router)
api_router.include_router(tasks_router)
api_router.include_router(results_router)


# 添加健康检查端点
@api_router.get("/health")
async def health_check():
    """API 健康检查"""
    return {"status": "healthy", "version": "1.0.0"}
