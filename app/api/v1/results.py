"""
结果查询 API 路由
"""
from typing import Optional, Union
from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse

from ...schemas.task import TaskStatus
from ...schemas.result import (
    OptimizationResult, StabilityResult,
    BulkModulusResult, HeatCapacityResult,
    QMOFEnergyResult, InteractionEnergyResult
)
from ...dependencies import get_task_service


router = APIRouter(prefix="/results", tags=["Results"])


ResultType = Union[
    OptimizationResult,
    StabilityResult,
    BulkModulusResult,
    HeatCapacityResult,
    QMOFEnergyResult,
    InteractionEnergyResult
]


@router.get("/{task_id}")
async def get_result(task_id: str):
    """
    获取任务结果
    
    返回结果的格式取决于任务类型
    """
    service = get_task_service()
    
    # 检查任务是否存在
    task = await service.get_task(task_id)
    if task is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Task '{task_id}' not found"
        )
    
    # 检查任务状态
    if task.status == TaskStatus.PENDING:
        return {
            "task_id": task_id,
            "status": "pending",
            "message": "Task is waiting to be executed"
        }
    
    if task.status == TaskStatus.RUNNING:
        return {
            "task_id": task_id,
            "status": "running",
            "progress": task.progress,
            "message": "Task is still running"
        }
    
    if task.status == TaskStatus.FAILED:
        return {
            "task_id": task_id,
            "status": "failed",
            "error": task.error_message
        }
    
    if task.status == TaskStatus.CANCELLED:
        return {
            "task_id": task_id,
            "status": "cancelled"
        }
    
    # 获取结果
    result = await service.get_result(task_id)
    if result is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Result for task '{task_id}' not found"
        )
    
    return result


@router.get("/{task_id}/download/{file_type}")
async def download_result_file(
    task_id: str,
    file_type: str
):
    """
    下载结果文件
    
    file_type: structure / trajectory / data
    """
    service = get_task_service()
    
    task = await service.get_task(task_id)
    if task is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Task '{task_id}' not found"
        )
    
    if task.status != TaskStatus.COMPLETED:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Task not completed"
        )
    
    result = await service.get_result(task_id)
    if result is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result not found"
        )
    
    # 根据文件类型确定路径
    file_path = None
    media_type = "application/octet-stream"
    
    if file_type == "structure":
        if hasattr(result, 'optimized_structure_file'):
            file_path = result.optimized_structure_file
            media_type = "chemical/x-cif"
        elif hasattr(result, 'final_structure_file'):
            file_path = result.final_structure_file
            media_type = "chemical/x-cif"
    
    elif file_type == "trajectory":
        if hasattr(result, 'trajectory_file'):
            file_path = result.trajectory_file
    
    elif file_type == "data":
        if hasattr(result, 'data_file'):
            file_path = result.data_file
            media_type = "application/json"
    
    if file_path is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"File type '{file_type}' not available for this result"
        )
    
    from pathlib import Path
    if not Path(file_path).exists():
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found on server"
        )
    
    return FileResponse(
        path=file_path,
        media_type=media_type,
        filename=Path(file_path).name
    )


@router.get("/{task_id}/summary")
async def get_result_summary(task_id: str):
    """获取结果摘要（适用于快速预览）"""
    service = get_task_service()
    
    task = await service.get_task(task_id)
    if task is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Task '{task_id}' not found"
        )
    
    if task.status != TaskStatus.COMPLETED:
        return {
            "task_id": task_id,
            "status": task.status.value,
            "completed": False
        }
    
    result = await service.get_result(task_id)
    if result is None:
        return {
            "task_id": task_id,
            "status": "completed",
            "result": None
        }
    
    # 根据结果类型提取关键信息
    summary = {
        "task_id": task_id,
        "task_type": task.task_type.value,
        "status": "completed"
    }
    
    if isinstance(result, OptimizationResult):
        summary.update({
            "final_energy": result.final_energy,
            "energy_per_atom": result.energy_per_atom,
            "converged": result.converged,
            "num_steps": result.num_steps
        })
    
    elif isinstance(result, StabilityResult):
        summary.update({
            "final_rmsd": result.final_rmsd,
            "volume_drift_percent": result.volume_drift_percent,
            "simulation_time_ps": result.simulation_time_ps
        })
    
    elif isinstance(result, BulkModulusResult):
        summary.update({
            "bulk_modulus": result.bulk_modulus,
            "r_squared": result.r_squared
        })
    
    elif isinstance(result, HeatCapacityResult):
        summary.update({
            "cv_300k": result.cv_300k,
            "has_imaginary_modes": result.has_imaginary_modes
        })
    
    return summary
