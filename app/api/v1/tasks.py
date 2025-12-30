"""
任务管理 API 路由
"""
from typing import List, Optional
from fastapi import APIRouter, HTTPException, status, BackgroundTasks

from ...schemas.task import (
    TaskType, TaskStatus, TaskInfo,
    OptimizationRequest, StabilityRequest,
    BulkModulusRequest, HeatCapacityRequest
)
from ...schemas.model import AVAILABLE_MODELS
from ...dependencies import get_task_service


router = APIRouter(prefix="/tasks", tags=["Tasks"])


def validate_model_key(model_key: str) -> None:
    """验证模型 key 是否有效"""
    if model_key not in AVAILABLE_MODELS:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid model key: '{model_key}'"
        )


# ========== 任务提交 ==========

@router.post("/optimization", response_model=TaskInfo)
async def submit_optimization_task(
    request: OptimizationRequest,
    background_tasks: BackgroundTasks
):
    """
    提交结构优化任务
    
    使用 BFGS/FIRE/LBFGS 优化器进行结构优化
    """
    validate_model_key(request.model_key)
    
    service = get_task_service()
    task = await service.create_optimization_task(request)
    
    # 在后台运行任务
    background_tasks.add_task(
        service.run_optimization,
        task.task_id,
        request
    )
    
    return task


@router.post("/stability", response_model=TaskInfo)
async def submit_stability_task(
    request: StabilityRequest,
    background_tasks: BackgroundTasks
):
    """
    提交 MD 稳定性模拟任务
    
    执行 NVT 平衡 + NPT 生产运行
    """
    validate_model_key(request.model_key)
    
    service = get_task_service()
    task = await service.create_stability_task(request)
    
    background_tasks.add_task(
        service.run_stability,
        task.task_id,
        request
    )
    
    return task


@router.post("/bulk-modulus", response_model=TaskInfo)
async def submit_bulk_modulus_task(
    request: BulkModulusRequest,
    background_tasks: BackgroundTasks
):
    """
    提交体积模量计算任务
    
    使用 E-V 曲线拟合 Birch-Murnaghan 状态方程
    """
    validate_model_key(request.model_key)
    
    service = get_task_service()
    task = await service.create_bulk_modulus_task(request)
    
    # 简化实现：创建任务并标记为 pending
    # 完整实现需要添加后台执行逻辑
    
    return task


@router.post("/heat-capacity", response_model=TaskInfo)
async def submit_heat_capacity_task(
    request: HeatCapacityRequest,
    background_tasks: BackgroundTasks
):
    """
    提交热容计算任务
    
    使用 Phonopy 进行声子计算
    """
    validate_model_key(request.model_key)
    
    service = get_task_service()
    task = await service.create_heat_capacity_task(request)
    
    return task


# ========== 任务查询 ==========

@router.get("/", response_model=List[TaskInfo])
async def list_tasks(
    task_type: Optional[TaskType] = None,
    status: Optional[TaskStatus] = None,
    model_key: Optional[str] = None,
    skip: int = 0,
    limit: int = 100
):
    """列出任务（支持筛选）"""
    service = get_task_service()
    tasks = await service.get_all_tasks()
    
    # 筛选
    if task_type:
        tasks = [t for t in tasks if t.task_type == task_type]
    if status:
        tasks = [t for t in tasks if t.status == status]
    if model_key:
        tasks = [t for t in tasks if t.model_key == model_key]
    
    return tasks[skip:skip + limit]


@router.get("/{task_id}", response_model=TaskInfo)
async def get_task(task_id: str):
    """获取任务状态和信息"""
    service = get_task_service()
    task = await service.get_task(task_id)
    
    if task is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Task '{task_id}' not found"
        )
    
    return task


@router.post("/{task_id}/cancel")
async def cancel_task(task_id: str):
    """取消任务"""
    service = get_task_service()
    
    success = await service.cancel_task(task_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot cancel task (not found or already completed)"
        )
    
    return {"status": "cancelled", "task_id": task_id}
