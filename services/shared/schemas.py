"""
MIRA 微服务共享 Schemas

所有模型 Worker 服务使用统一的请求/响应格式
"""
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
from enum import Enum


# ============================================
# 枚举类型
# ============================================
class TaskType(str, Enum):
    """任务类型"""
    SINGLE_POINT = "single_point"
    OPTIMIZATION = "optimization"
    STABILITY = "stability"
    BULK_MODULUS = "bulk_modulus"
    HEAT_CAPACITY = "heat_capacity"
    MD_SIMULATION = "md_simulation"


class TaskStatus(str, Enum):
    """任务状态"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


# ============================================
# 原子结构
# ============================================
class AtomData(BaseModel):
    """原子数据"""
    symbols: List[str] = Field(..., description="元素符号列表")
    positions: List[List[float]] = Field(..., description="原子位置 (Å)")
    cell: Optional[List[List[float]]] = Field(None, description="晶胞矩阵 3x3")
    pbc: List[bool] = Field(default=[True, True, True], description="周期性边界")


# ============================================
# 单点计算
# ============================================
class SinglePointRequest(BaseModel):
    """单点能量计算请求"""
    atoms: AtomData
    model_name: str = Field(..., description="模型名称")
    compute_stress: bool = Field(default=True, description="是否计算应力")
    compute_forces: bool = Field(default=True, description="是否计算力")


class SinglePointResponse(BaseModel):
    """单点能量计算响应"""
    energy: float = Field(..., description="总能量 (eV)")
    forces: Optional[List[List[float]]] = Field(None, description="原子力 (eV/Å)")
    stress: Optional[List[float]] = Field(None, description="应力张量 (eV/Å³)")
    energy_per_atom: float = Field(..., description="每原子能量 (eV/atom)")
    max_force: Optional[float] = Field(None, description="最大原子力 (eV/Å)")


# ============================================
# 结构优化
# ============================================
class OptimizationRequest(BaseModel):
    """结构优化请求"""
    atoms: AtomData
    model_name: str = Field(..., description="模型名称")
    fmax: float = Field(default=0.05, description="力收敛阈值 (eV/Å)")
    max_steps: int = Field(default=500, description="最大优化步数")
    optimizer: str = Field(default="BFGS", description="优化器类型")
    use_d3: bool = Field(default=True, description="是否使用 D3 色散校正")
    fix_cell: bool = Field(default=False, description="是否固定晶胞")


class OptimizationResponse(BaseModel):
    """结构优化响应"""
    converged: bool = Field(..., description="是否收敛")
    steps: int = Field(..., description="优化步数")
    initial_energy: float = Field(..., description="初始能量 (eV)")
    final_energy: float = Field(..., description="最终能量 (eV)")
    energy_change: float = Field(..., description="能量变化 (eV)")
    final_forces_max: float = Field(..., description="最终最大力 (eV/Å)")
    optimized_atoms: AtomData = Field(..., description="优化后的结构")


# ============================================
# 稳定性分析
# ============================================
class StabilityRequest(BaseModel):
    """稳定性分析请求"""
    atoms: AtomData
    model_name: str = Field(..., description="模型名称")
    temperature: float = Field(default=300.0, description="温度 (K)")
    pressure: float = Field(default=0.0, description="压力 (bar)")
    timestep: float = Field(default=1.0, description="时间步长 (fs)")
    equilibration_steps: int = Field(default=1000, description="平衡步数")
    production_steps: int = Field(default=5000, description="采样步数")
    use_d3: bool = Field(default=True, description="是否使用 D3 色散校正")


class StabilityResponse(BaseModel):
    """稳定性分析响应"""
    is_stable: bool = Field(..., description="是否稳定")
    mean_energy: float = Field(..., description="平均能量 (eV)")
    energy_std: float = Field(..., description="能量标准差 (eV)")
    mean_temperature: float = Field(..., description="平均温度 (K)")
    mean_pressure: float = Field(..., description="平均压力 (bar)")
    volume_change_percent: float = Field(..., description="体积变化百分比")
    max_displacement: float = Field(..., description="最大原子位移 (Å)")
    trajectory_length: int = Field(..., description="轨迹长度")


# ============================================
# 体积模量
# ============================================
class BulkModulusRequest(BaseModel):
    """体积模量计算请求"""
    atoms: AtomData
    model_name: str = Field(..., description="模型名称")
    strain_range: float = Field(default=0.06, description="应变范围")
    num_points: int = Field(default=7, description="采样点数")
    use_d3: bool = Field(default=True, description="是否使用 D3 色散校正")


class BulkModulusResponse(BaseModel):
    """体积模量计算响应"""
    bulk_modulus: float = Field(..., description="体积模量 (GPa)")
    equilibrium_volume: float = Field(..., description="平衡体积 (Å³)")
    equilibrium_energy: float = Field(..., description="平衡能量 (eV)")
    fitting_r2: float = Field(..., description="拟合 R² 值")
    volumes: List[float] = Field(..., description="体积数据点")
    energies: List[float] = Field(..., description="能量数据点")


# ============================================
# 热容
# ============================================
class HeatCapacityRequest(BaseModel):
    """热容计算请求"""
    atoms: AtomData
    model_name: str = Field(..., description="模型名称")
    temperatures: List[float] = Field(
        default=[100, 200, 300, 400, 500],
        description="温度列表 (K)"
    )
    supercell: List[int] = Field(default=[2, 2, 2], description="超胞大小")
    use_d3: bool = Field(default=True, description="是否使用 D3 色散校正")


class HeatCapacityResponse(BaseModel):
    """热容计算响应"""
    temperatures: List[float] = Field(..., description="温度 (K)")
    heat_capacities: List[float] = Field(..., description="热容 Cv (J/mol·K)")
    free_energies: List[float] = Field(..., description="自由能 (eV)")
    entropies: List[float] = Field(..., description="熵 (J/mol·K)")


# ============================================
# 通用任务请求/响应
# ============================================
class TaskRequest(BaseModel):
    """通用任务请求"""
    task_id: str = Field(..., description="任务 ID")
    task_type: TaskType = Field(..., description="任务类型")
    payload: Dict[str, Any] = Field(..., description="任务参数")


class TaskResponse(BaseModel):
    """通用任务响应"""
    task_id: str = Field(..., description="任务 ID")
    status: TaskStatus = Field(..., description="任务状态")
    result: Optional[Dict[str, Any]] = Field(None, description="任务结果")
    error: Optional[str] = Field(None, description="错误信息")
    execution_time: Optional[float] = Field(None, description="执行时间 (秒)")


# ============================================
# 健康检查
# ============================================
class HealthResponse(BaseModel):
    """健康检查响应"""
    status: str = Field(default="healthy")
    service: str = Field(..., description="服务名称")
    models: List[str] = Field(..., description="可用模型列表")
    device: str = Field(..., description="计算设备")
    gpu_available: bool = Field(..., description="GPU 是否可用")
    gpu_name: Optional[str] = Field(None, description="GPU 名称")
    memory_usage_mb: Optional[float] = Field(None, description="显存使用 (MB)")


class ModelsListResponse(BaseModel):
    """模型列表响应"""
    models: List[str] = Field(..., description="可用模型列表")
    service: str = Field(..., description="服务名称")
