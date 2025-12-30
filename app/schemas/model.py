"""
模型相关的 Pydantic 数据模型
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum


class ModelFamily(str, Enum):
    """模型系列枚举"""
    MACE = "mace"
    ORB = "orb"
    OMAT24 = "omat24"
    GRACE = "grace"
    MATTERSIM = "mattersim"
    SEVENNET = "sevennet"
    POSEGNN = "posegnn"
    MATGL = "matgl"


class ModelInfo(BaseModel):
    """模型信息"""
    key: str = Field(..., description="模型唯一标识符，如 'mace_prod_b3'")
    name: str = Field(..., description="模型显示名称，如 'MACE-MP-0b3-medium'")
    family: ModelFamily = Field(..., description="模型系列")
    description: Optional[str] = Field(None, description="模型描述")
    version: Optional[str] = Field(None, description="模型版本")
    supports_d3: bool = Field(True, description="是否支持 D3 色散校正")
    is_loaded: bool = Field(False, description="是否已加载到内存")
    memory_usage_mb: Optional[float] = Field(None, description="内存占用 (MB)")
    device: Optional[str] = Field(None, description="当前运行设备")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "key": "mace_prod_b3",
                "name": "MACE-MP-0b3-medium",
                "family": "mace",
                "description": "MACE medium model trained on Materials Project",
                "supports_d3": True,
                "is_loaded": False
            }
        }
    }


class ModelLoadRequest(BaseModel):
    """模型加载请求"""
    device: str = Field("cuda", description="运行设备: cuda / cpu")
    dtype: str = Field("float32", description="数据类型: float32 / float64")
    enable_d3: bool = Field(True, description="是否启用 D3 色散校正")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "device": "cuda",
                "dtype": "float32",
                "enable_d3": True
            }
        }
    }


class ModelLoadResponse(BaseModel):
    """模型加载响应"""
    success: bool
    model_key: str
    message: str
    memory_usage_mb: Optional[float] = None


class ModelStatusResponse(BaseModel):
    """模型状态响应"""
    loaded_models: List[str]
    total_memory_mb: float
    available_models: List[str]


# 预定义的模型配置
AVAILABLE_MODELS: Dict[str, Dict[str, Any]] = {
    # MACE 系列
    "mace_prod_b3": {
        "name": "MACE-MP-0b3-medium",
        "family": ModelFamily.MACE,
        "description": "MACE medium model (MP-0b3)",
        "supports_d3": True
    },
    "mace_prod": {
        "name": "MACE-MPA-0-medium",
        "family": ModelFamily.MACE,
        "description": "MACE MPA medium model",
        "supports_d3": True
    },
    "mace_prod_omat": {
        "name": "MACE-OMAT-0-medium",
        "family": ModelFamily.MACE,
        "description": "MACE OMAT medium model",
        "supports_d3": True
    },
    "mace_prod_matpes": {
        "name": "MACE-MatPES-r2scan",
        "family": ModelFamily.MACE,
        "description": "MACE MatPES r2SCAN model",
        "supports_d3": True
    },
    "mace_prod_mof": {
        "name": "MACE-MOF-v1",
        "family": ModelFamily.MACE,
        "description": "MACE model trained for MOF",
        "supports_d3": True
    },
    
    # ORB 系列
    "orb_prod_mp": {
        "name": "orb-mptraj-only-v2",
        "family": ModelFamily.ORB,
        "description": "ORB model trained on MP trajectory",
        "supports_d3": False
    },
    "orb_prod": {
        "name": "orb-d3-v2",
        "family": ModelFamily.ORB,
        "description": "ORB model with D3 correction",
        "supports_d3": True
    },
    "orb3": {
        "name": "orb-v3-direct-20-omat",
        "family": ModelFamily.ORB,
        "description": "ORB v3 direct model",
        "supports_d3": True
    },
    "orb_prod_v3": {
        "name": "orb-v3-conservative-inf-omat",
        "family": ModelFamily.ORB,
        "description": "ORB v3 conservative OMAT model",
        "supports_d3": True
    },
    "orb_prod_v3_mp": {
        "name": "orb-v3-conservative-inf-mpa",
        "family": ModelFamily.ORB,
        "description": "ORB v3 conservative MPA model",
        "supports_d3": True
    },
    
    # OMAT24 (FAIRChem) 系列
    "omat24_prod_mp": {
        "name": "eqV2_dens_86M_mp",
        "family": ModelFamily.OMAT24,
        "description": "EquiformerV2 86M MP model",
        "supports_d3": True
    },
    "omat24_prod": {
        "name": "eqV2_86M_omat_mp_salex",
        "family": ModelFamily.OMAT24,
        "description": "EquiformerV2 86M OMAT model",
        "supports_d3": True
    },
    "omat24_prod_esen": {
        "name": "esen_30m_oam",
        "family": ModelFamily.OMAT24,
        "description": "eSEN 30M OAM model",
        "supports_d3": True
    },
    "omat24_prod_esen_mp": {
        "name": "esen_30m_mptrj",
        "family": ModelFamily.OMAT24,
        "description": "eSEN 30M MP trajectory model",
        "supports_d3": True
    },
    
    # GRACE 系列
    "grace_prod": {
        "name": "GRACE-2L-MP-r6",
        "family": ModelFamily.GRACE,
        "description": "GRACE 2-layer MP model",
        "supports_d3": True
    },
    "grace_prod_oam": {
        "name": "GRACE-2L-OAM",
        "family": ModelFamily.GRACE,
        "description": "GRACE 2-layer OAM model",
        "supports_d3": True
    },
    "grace_prod_omat": {
        "name": "GRACE-2L-OMAT",
        "family": ModelFamily.GRACE,
        "description": "GRACE 2-layer OMAT model",
        "supports_d3": True
    },
    
    # MatterSim
    "mattersim_prod": {
        "name": "MatterSim-v1.0.0-5M",
        "family": ModelFamily.MATTERSIM,
        "description": "MatterSim 5M model",
        "supports_d3": True
    },
    
    # SevenNet
    "sevennet_prod": {
        "name": "7net-0",
        "family": ModelFamily.SEVENNET,
        "description": "SevenNet base model",
        "supports_d3": True
    },
    "sevennet_prod_l3i5": {
        "name": "7net-l3i5",
        "family": ModelFamily.SEVENNET,
        "description": "SevenNet L3I5 model",
        "supports_d3": True
    },
    "sevennet_prod_ompa": {
        "name": "7net-mf-ompa",
        "family": ModelFamily.SEVENNET,
        "description": "SevenNet MF OMPA model",
        "supports_d3": True
    },
    
    # PosEGNN (IBM)
    "posegnn_prod": {
        "name": "pos-egnn.v1-6M",
        "family": ModelFamily.POSEGNN,
        "description": "IBM PosEGNN 6M model",
        "supports_d3": True
    },
}
