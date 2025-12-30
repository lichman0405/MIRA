"""
结构相关的 Pydantic 数据模型
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime
import uuid


class StructureUpload(BaseModel):
    """结构上传请求"""
    name: str = Field(..., description="结构名称")
    format: str = Field(..., description="文件格式: cif / poscar / xyz / json")
    content: str = Field(..., description="文件内容或 Base64 编码")
    metadata: Optional[Dict[str, Any]] = Field(None, description="附加元数据")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "name": "MOF-5",
                "format": "cif",
                "content": "data_MOF-5\n_cell_length_a   25.832...",
                "metadata": {"source": "CSD", "refcode": "SAHYOG"}
            }
        }
    }


class StructureInfo(BaseModel):
    """结构信息"""
    id: str = Field(default_factory=lambda: str(uuid.uuid4()), description="结构唯一标识符")
    name: str = Field(..., description="结构名称")
    formula: str = Field(..., description="化学式")
    num_atoms: int = Field(..., description="原子数量")
    cell_volume: float = Field(..., description="晶胞体积 (Å³)")
    cell_parameters: Dict[str, float] = Field(..., description="晶胞参数 a, b, c, α, β, γ")
    space_group: Optional[str] = Field(None, description="空间群")
    elements: List[str] = Field(..., description="包含的元素")
    file_path: str = Field(..., description="存储路径")
    created_at: datetime = Field(default_factory=datetime.now, description="创建时间")
    metadata: Optional[Dict[str, Any]] = Field(None, description="附加元数据")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "id": "550e8400-e29b-41d4-a716-446655440000",
                "name": "MOF-5",
                "formula": "Zn4O(BDC)3",
                "num_atoms": 424,
                "cell_volume": 17237.5,
                "cell_parameters": {
                    "a": 25.832, "b": 25.832, "c": 25.832,
                    "alpha": 90.0, "beta": 90.0, "gamma": 90.0
                },
                "space_group": "Fm-3m",
                "elements": ["Zn", "O", "C", "H"],
                "file_path": "./data/structures/550e8400.cif",
                "created_at": "2024-01-15T10:30:00"
            }
        }
    }


class StructureListResponse(BaseModel):
    """结构列表响应"""
    total: int = Field(..., description="总数")
    structures: List[StructureInfo] = Field(..., description="结构列表")


class StructurePreview(BaseModel):
    """结构预览数据（用于可视化）"""
    id: str
    name: str
    formula: str
    positions: List[List[float]] = Field(..., description="原子坐标 [[x, y, z], ...]")
    symbols: List[str] = Field(..., description="原子符号")
    cell: List[List[float]] = Field(..., description="晶胞矩阵 3x3")
    pbc: List[bool] = Field(..., description="周期性边界条件 [x, y, z]")
