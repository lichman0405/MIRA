"""
任务管理服务
"""
from typing import Dict, Optional, List, Any
from datetime import datetime
from pathlib import Path
import uuid
import asyncio
from concurrent.futures import ThreadPoolExecutor

from ase import Atoms

from ..schemas.task import (
    TaskType, TaskStatus, TaskInfo,
    OptimizationRequest, StabilityRequest,
    BulkModulusRequest, HeatCapacityRequest
)
from ..schemas.result import (
    OptimizationResult, StabilityResult,
    BulkModulusResult, HeatCapacityResult
)
from ..models.registry import ModelRegistry
from ..core.ase_utils import load_structure
from .optimization import OptimizationService
from .stability import StabilityService
from .bulk_modulus import BulkModulusService
from .heat_capacity import HeatCapacityService


class TaskService:
    """任务管理服务"""
    
    def __init__(
        self,
        model_registry: ModelRegistry,
        structures_dir: Path,
        results_dir: Path,
        max_workers: int = 4
    ):
        self.model_registry = model_registry
        self.structures_dir = Path(structures_dir)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # 任务存储（生产环境应使用数据库）
        self._tasks: Dict[str, TaskInfo] = {}
        self._results: Dict[str, Any] = {}
        
        # 线程池用于运行计算任务
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        
        # 初始化服务
        self._opt_service = OptimizationService(results_dir)
        self._stability_service = StabilityService(results_dir)
        self._bulk_service = BulkModulusService(results_dir)
        self._heat_service = HeatCapacityService(results_dir)
    
    def _load_structure(self, structure_id: str) -> Atoms:
        """加载结构文件"""
        # 尝试常见格式
        for ext in ['.cif', '.vasp', '.xyz', '.json']:
            path = self.structures_dir / f"{structure_id}{ext}"
            if path.exists():
                return load_structure(path)
        
        # 尝试直接使用 structure_id 作为文件名
        path = self.structures_dir / structure_id
        if path.exists():
            return load_structure(path)
        
        raise FileNotFoundError(f"Structure not found: {structure_id}")
    
    def _update_task(
        self,
        task_id: str,
        status: Optional[TaskStatus] = None,
        progress: Optional[float] = None,
        error_message: Optional[str] = None,
        result_path: Optional[str] = None
    ) -> None:
        """更新任务状态"""
        if task_id not in self._tasks:
            return
        
        task = self._tasks[task_id]
        
        if status is not None:
            task.status = status
            if status == TaskStatus.RUNNING and task.started_at is None:
                task.started_at = datetime.now()
            elif status in [TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED]:
                task.completed_at = datetime.now()
        
        if progress is not None:
            task.progress = progress
        
        if error_message is not None:
            task.error_message = error_message
        
        if result_path is not None:
            task.result_path = result_path
    
    # ========== 任务创建 ==========
    
    async def create_optimization_task(
        self,
        request: OptimizationRequest
    ) -> TaskInfo:
        """创建优化任务"""
        task = TaskInfo(
            task_id=str(uuid.uuid4()),
            task_type=TaskType.OPTIMIZATION,
            status=TaskStatus.PENDING,
            model_key=request.model_key,
            structure_id=request.structure_id
        )
        self._tasks[task.task_id] = task
        return task
    
    async def create_stability_task(
        self,
        request: StabilityRequest
    ) -> TaskInfo:
        """创建稳定性模拟任务"""
        task = TaskInfo(
            task_id=str(uuid.uuid4()),
            task_type=TaskType.STABILITY,
            status=TaskStatus.PENDING,
            model_key=request.model_key,
            structure_id=request.structure_id
        )
        self._tasks[task.task_id] = task
        return task
    
    async def create_bulk_modulus_task(
        self,
        request: BulkModulusRequest
    ) -> TaskInfo:
        """创建体积模量计算任务"""
        task = TaskInfo(
            task_id=str(uuid.uuid4()),
            task_type=TaskType.BULK_MODULUS,
            status=TaskStatus.PENDING,
            model_key=request.model_key,
            structure_id=request.structure_id
        )
        self._tasks[task.task_id] = task
        return task
    
    async def create_heat_capacity_task(
        self,
        request: HeatCapacityRequest
    ) -> TaskInfo:
        """创建热容计算任务"""
        task = TaskInfo(
            task_id=str(uuid.uuid4()),
            task_type=TaskType.HEAT_CAPACITY,
            status=TaskStatus.PENDING,
            model_key=request.model_key,
            structure_id=request.structure_id
        )
        self._tasks[task.task_id] = task
        return task
    
    # ========== 任务执行 ==========
    
    def _run_optimization_sync(
        self,
        task_id: str,
        request: OptimizationRequest
    ) -> OptimizationResult:
        """同步执行优化任务"""
        self._update_task(task_id, status=TaskStatus.RUNNING)
        
        try:
            # 加载结构和模型
            atoms = self._load_structure(request.structure_id)
            model = self.model_registry.get_or_load_model(request.model_key)
            
            # 进度回调
            def progress_cb(p):
                self._update_task(task_id, progress=p)
            
            # 执行优化
            result = self._opt_service.run(
                atoms=atoms,
                model=model,
                request=request,
                task_id=task_id,
                progress_callback=progress_cb
            )
            
            self._results[task_id] = result
            self._update_task(
                task_id,
                status=TaskStatus.COMPLETED,
                progress=100.0,
                result_path=result.optimized_structure_file
            )
            
            return result
            
        except Exception as e:
            self._update_task(
                task_id,
                status=TaskStatus.FAILED,
                error_message=str(e)
            )
            raise
    
    async def run_optimization(
        self,
        task_id: str,
        request: Optional[OptimizationRequest] = None
    ) -> OptimizationResult:
        """异步执行优化任务"""
        if request is None:
            # 从存储的任务中获取请求（简化实现）
            raise ValueError("Request is required")
        
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            self._executor,
            self._run_optimization_sync,
            task_id,
            request
        )
    
    def _run_stability_sync(
        self,
        task_id: str,
        request: StabilityRequest
    ) -> StabilityResult:
        """同步执行稳定性模拟"""
        self._update_task(task_id, status=TaskStatus.RUNNING)
        
        try:
            atoms = self._load_structure(request.structure_id)
            model = self.model_registry.get_or_load_model(request.model_key)
            
            def progress_cb(p):
                self._update_task(task_id, progress=p)
            
            result = self._stability_service.run(
                atoms=atoms,
                model=model,
                request=request,
                task_id=task_id,
                progress_callback=progress_cb
            )
            
            self._results[task_id] = result
            self._update_task(
                task_id,
                status=TaskStatus.COMPLETED,
                progress=100.0,
                result_path=result.trajectory_file
            )
            
            return result
            
        except Exception as e:
            self._update_task(
                task_id,
                status=TaskStatus.FAILED,
                error_message=str(e)
            )
            raise
    
    async def run_stability(
        self,
        task_id: str,
        request: StabilityRequest
    ) -> StabilityResult:
        """异步执行稳定性模拟"""
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            self._executor,
            self._run_stability_sync,
            task_id,
            request
        )
    
    # ========== 任务查询 ==========
    
    async def get_task(self, task_id: str) -> Optional[TaskInfo]:
        """获取任务信息"""
        return self._tasks.get(task_id)
    
    async def get_all_tasks(self) -> List[TaskInfo]:
        """获取所有任务"""
        return list(self._tasks.values())
    
    async def get_result(self, task_id: str) -> Optional[Any]:
        """获取任务结果"""
        return self._results.get(task_id)
    
    async def cancel_task(self, task_id: str) -> bool:
        """取消任务"""
        if task_id not in self._tasks:
            return False
        
        task = self._tasks[task_id]
        if task.status not in [TaskStatus.PENDING, TaskStatus.RUNNING]:
            return False
        
        self._update_task(task_id, status=TaskStatus.CANCELLED)
        return True
