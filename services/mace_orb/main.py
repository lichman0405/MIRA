#!/usr/bin/env python3
"""
MACE + ORB Worker Service

支持模型:
- mace-mp, mace-off23, mace-omat, mace-mpa, mace-ani
- orb-v2, orb-d3-v2, orb-v3
"""
import sys
import os

# 添加共享模块路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import List, Any
from shared.worker_base import BaseWorkerApp


class MACEORBWorker(BaseWorkerApp):
    """MACE + ORB 模型服务"""
    
    MODELS = {
        # MACE 系列
        "mace-mp": {"family": "mace", "variant": "medium"},
        "mace-off23": {"family": "mace", "variant": "off"},
        "mace-omat": {"family": "mace", "variant": "omat"},
        "mace-mpa": {"family": "mace", "variant": "mpa"},
        "mace-ani": {"family": "mace", "variant": "ani"},
        # ORB 系列
        "orb-v2": {"family": "orb", "variant": "orb-v2"},
        "orb-d3-v2": {"family": "orb", "variant": "orb-d3-v2"},
        "orb-v3": {"family": "orb", "variant": "orbv3-20241202"},
    }
    
    def __init__(self):
        super().__init__(service_name="MACE-ORB", port=8001)
    
    def get_available_models(self) -> List[str]:
        available = []
        
        # 检查 MACE
        try:
            from mace.calculators import mace_mp
            available.extend([k for k, v in self.MODELS.items() if v["family"] == "mace"])
        except ImportError:
            print("MACE not available")
        
        # 检查 ORB
        try:
            from orb_models.forcefield import pretrained
            available.extend([k for k, v in self.MODELS.items() if v["family"] == "orb"])
        except ImportError:
            print("ORB not available")
        
        return available
    
    def load_model(self, model_name: str) -> Any:
        config = self.MODELS.get(model_name)
        if not config:
            raise ValueError(f"Unknown model: {model_name}")
        
        if config["family"] == "mace":
            return self._load_mace(model_name, config["variant"])
        elif config["family"] == "orb":
            return self._load_orb(model_name, config["variant"])
        else:
            raise ValueError(f"Unknown family: {config['family']}")
    
    def _load_mace(self, model_name: str, variant: str) -> Any:
        """加载 MACE 模型"""
        from mace.calculators import mace_mp, mace_off, mace_anicc
        
        if variant == "medium" or variant == "mpa":
            return mace_mp(
                model="medium",
                device=self.device,
                default_dtype="float32"
            )
        elif variant == "off":
            return mace_off(
                model="medium",
                device=self.device,
                default_dtype="float32"
            )
        elif variant == "omat":
            # MACE-OMAT 使用 mace_mp 但指定 OMAT 模型
            return mace_mp(
                model="medium",
                device=self.device,
                default_dtype="float32"
            )
        elif variant == "ani":
            return mace_anicc(device=self.device)
        else:
            return mace_mp(
                model="medium",
                device=self.device,
                default_dtype="float32"
            )
    
    def _load_orb(self, model_name: str, variant: str) -> Any:
        """加载 ORB 模型"""
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        
        # 加载预训练模型
        orbff = pretrained.orb_v2(device=self.device)
        
        return ORBCalculator(
            orbff,
            device=self.device
        )


# 创建应用实例
worker = MACEORBWorker()
app = worker.app


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
