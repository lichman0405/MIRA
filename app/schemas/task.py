"""
任务相关的 Pydantic 数据模型
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum
from datetime import datetime
import uuid


class TaskType(str, Enum):
    """任务类型枚举"""
    OPTIMIZATION = "optimization"
    STABILITY = "stability"
    BULK_MODULUS = "bulk_modulus"
    HEAT_CAPACITY = "heat_capacity"
    QMOF_ENERGY = "qmof_energy"
    INTERACTION_ENERGY = "interaction_energy"


class TaskStatus(str, Enum):
    """任务状态枚举"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class TaskInfo(BaseModel):
    """任务信息"""
    task_id: str = Field(default_factory=lambda: str(uuid.uuid4()), description="任务唯一标识符")
    task_type: TaskType = Field(..., description="任务类型")
    status: TaskStatus = Field(default=TaskStatus.PENDING, description="任务状态")
    progress: float = Field(default=0.0, ge=0.0, le=100.0, description="任务进度 (0-100)")
    model_key: str = Field(..., description="使用的模型")
    structure_id: str = Field(..., description="结构 ID")
    created_at: datetime = Field(default_factory=datetime.now, description="创建时间")
    started_at: Optional[datetime] = Field(None, description="开始时间")
    completed_at: Optional[datetime] = Field(None, description="完成时间")
    error_message: Optional[str] = Field(None, description="错误信息")
    result_path: Optional[str] = Field(None, description="结果文件路径")


# ========== 结构优化任务 ==========
class OptimizationRequest(BaseModel):
    """结构优化任务请求"""
    structure_id: str = Field(..., description="结构 ID")
    model_key: str = Field(..., description="使用的模型")
    fmax: float = Field(0.01, gt=0, description="力收敛标准 (eV/Å)")
    max_steps: int = Field(500, gt=0, le=10000, description="最大优化步数")
    optimizer: str = Field("BFGS", description="优化器: BFGS / FIRE / LBFGS")
    use_filter: bool = Field(True, description="使用 FrechetCellFilter 进行全松弛")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_id": "550e8400-e29b-41d4-a716-446655440000",
                "model_key": "mace_prod_b3",
                "fmax": 0.01,
                "max_steps": 500,
                "optimizer": "BFGS",
                "use_filter": True,
                "enable_d3": True
            }
        }
    }


# ========== MD 稳定性模拟 ==========
class StabilityRequest(BaseModel):
    """MD 稳定性模拟任务请求"""
    structure_id: str = Field(..., description="结构 ID")
    model_key: str = Field(..., description="使用的模型")
    
    # 预优化参数
    pre_optimize: bool = Field(True, description="是否先进行结构优化")
    fmax: float = Field(0.05, gt=0, description="预优化力收敛标准 (eV/Å)")
    
    # NVT 平衡参数
    nvt_temperature: float = Field(300.0, gt=0, description="NVT 温度 (K)")
    nvt_steps: int = Field(1000, gt=0, description="NVT 步数")
    
    # NPT 生产运行参数
    npt_temperature: float = Field(300.0, gt=0, description="NPT 温度 (K)")
    npt_pressure: float = Field(1.0, description="NPT 压力 (bar)")
    npt_steps: int = Field(10000, gt=0, description="NPT 步数")
    
    # 通用 MD 参数
    timestep: float = Field(1.0, gt=0, description="时间步长 (fs)")
    thermostat: str = Field("langevin", description="热浴类型: langevin / berendsen / nose-hoover")
    barostat: str = Field("berendsen", description="压力浴类型: berendsen / mtk")
    trajectory_interval: int = Field(10, gt=0, description="轨迹保存间隔")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_id": "550e8400-e29b-41d4-a716-446655440000",
                "model_key": "mace_prod_b3",
                "pre_optimize": True,
                "nvt_temperature": 300.0,
                "nvt_steps": 1000,
                "npt_temperature": 300.0,
                "npt_pressure": 1.0,
                "npt_steps": 10000,
                "timestep": 1.0,
                "thermostat": "langevin",
                "barostat": "berendsen",
                "enable_d3": True
            }
        }
    }


# ========== 体积模量计算 ==========
class BulkModulusRequest(BaseModel):
    """体积模量计算任务请求"""
    structure_id: str = Field(..., description="结构 ID")
    model_key: str = Field(..., description="使用的模型")
    volume_range: float = Field(0.1, gt=0, le=0.5, description="体积变化范围 (±比例)")
    num_points: int = Field(11, ge=5, le=51, description="采样点数")
    eos_type: str = Field("birch-murnaghan", description="状态方程类型")
    pre_optimize: bool = Field(True, description="是否先优化结构")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_id": "550e8400-e29b-41d4-a716-446655440000",
                "model_key": "mace_prod_b3",
                "volume_range": 0.1,
                "num_points": 11,
                "eos_type": "birch-murnaghan",
                "enable_d3": True
            }
        }
    }


# ========== 热容计算 ==========
class HeatCapacityRequest(BaseModel):
    """热容计算任务请求"""
    structure_id: str = Field(..., description="结构 ID")
    model_key: str = Field(..., description="使用的模型")
    supercell: List[int] = Field([2, 2, 2], description="超胞大小 [nx, ny, nz]")
    temperature_range: List[float] = Field([0, 1000], description="温度范围 [min, max] (K)")
    temperature_points: int = Field(101, gt=0, description="温度采样点数")
    mesh: List[int] = Field([8, 8, 8], description="Phonon k 点网格")
    displacement: float = Field(0.01, gt=0, description="有限位移大小 (Å)")
    pre_optimize: bool = Field(True, description="是否先优化结构")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_id": "550e8400-e29b-41d4-a716-446655440000",
                "model_key": "mace_prod_b3",
                "supercell": [2, 2, 2],
                "temperature_range": [0, 1000],
                "mesh": [8, 8, 8],
                "enable_d3": True
            }
        }
    }


# ========== QMOF 能量对比 ==========
class QMOFEnergyRequest(BaseModel):
    """QMOF 能量对比任务请求"""
    structure_ids: List[str] = Field(..., description="结构 ID 列表")
    model_key: str = Field(..., description="使用的模型")
    reference_energies: Optional[Dict[str, float]] = Field(None, description="DFT 参考能量")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_ids": ["id1", "id2", "id3"],
                "model_key": "mace_prod_b3",
                "reference_energies": {"id1": -100.5, "id2": -200.3},
                "enable_d3": True
            }
        }
    }


# ========== 相互作用能 ==========
class InteractionEnergyRequest(BaseModel):
    """相互作用能计算任务请求"""
    host_structure_id: str = Field(..., description="主体 MOF 结构 ID")
    guest_molecule: str = Field(..., description="客体分子: CO2 / H2O / N2 / CH4")
    model_key: str = Field(..., description="使用的模型")
    binding_sites: Optional[List[Dict[str, Any]]] = Field(None, description="吸附位点")
    num_configurations: int = Field(10, gt=0, description="采样构型数")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "host_structure_id": "550e8400-e29b-41d4-a716-446655440000",
                "guest_molecule": "CO2",
                "model_key": "mace_prod_b3",
                "num_configurations": 10,
                "enable_d3": True
            }
        }
    }


# ========== 批量基准测试 ==========
class BatchBenchmarkRequest(BaseModel):
    """批量基准测试请求"""
    structure_ids: List[str] = Field(..., description="结构 ID 列表")
    model_keys: List[str] = Field(..., description="模型列表")
    task_types: List[TaskType] = Field(..., description="任务类型列表")
    enable_d3: bool = Field(True, description="启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "structure_ids": ["id1", "id2"],
                "model_keys": ["mace_prod_b3", "orb_prod"],
                "task_types": ["optimization", "stability"],
                "enable_d3": True
            }
        }
    }


class BatchInfo(BaseModel):
    """批量任务信息"""
    batch_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    total_tasks: int
    completed_tasks: int = 0
    failed_tasks: int = 0
    task_ids: List[str] = Field(default_factory=list)
    status: TaskStatus = TaskStatus.PENDING
    created_at: datetime = Field(default_factory=datetime.now)
