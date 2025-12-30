#!/usr/bin/env python3
"""
FAIRChem + SevenNet Worker Service

支持模型:
- omat24-base, omat24-large, eqv2-omat, eqv2-mptrj
- sevennet-0, sevennet-mf-ompa, sevennet-l3i5
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import List, Any
from shared.worker_base import BaseWorkerApp


class FAIRChemSevenNetWorker(BaseWorkerApp):
    """FAIRChem + SevenNet 模型服务"""
    
    MODELS = {
        # FAIRChem/OMAT24 系列
        "omat24-base": {"family": "fairchem", "variant": "base"},
        "omat24-large": {"family": "fairchem", "variant": "large"},
        "eqv2-omat": {"family": "fairchem", "variant": "eqv2-omat"},
        "eqv2-mptrj": {"family": "fairchem", "variant": "eqv2-mptrj"},
        # SevenNet 系列
        "sevennet-0": {"family": "sevennet", "variant": "0"},
        "sevennet-mf-ompa": {"family": "sevennet", "variant": "mf-ompa"},
        "sevennet-l3i5": {"family": "sevennet", "variant": "l3i5"},
    }
    
    def __init__(self):
        super().__init__(service_name="FAIRChem-SevenNet", port=8002)
    
    def get_available_models(self) -> List[str]:
        available = []
        
        # 检查 FAIRChem
        try:
            from fairchem.core import OCPCalculator
            available.extend([k for k, v in self.MODELS.items() if v["family"] == "fairchem"])
        except ImportError:
            print("FAIRChem not available")
        
        # 检查 SevenNet
        try:
            from sevenn.sevennet_calculator import SevenNetCalculator
            available.extend([k for k, v in self.MODELS.items() if v["family"] == "sevennet"])
        except ImportError:
            print("SevenNet not available")
        
        return available
    
    def load_model(self, model_name: str) -> Any:
        config = self.MODELS.get(model_name)
        if not config:
            raise ValueError(f"Unknown model: {model_name}")
        
        if config["family"] == "fairchem":
            return self._load_fairchem(model_name, config["variant"])
        elif config["family"] == "sevennet":
            return self._load_sevennet(model_name, config["variant"])
        else:
            raise ValueError(f"Unknown family: {config['family']}")
    
    def _load_fairchem(self, model_name: str, variant: str) -> Any:
        """加载 FAIRChem 模型"""
        from fairchem.core import OCPCalculator
        
        # 模型路径映射
        model_paths = {
            "base": "EquiformerV2-31M-S2EF-OC20-All+MD",
            "large": "EquiformerV2-153M-S2EF-OC20-All+MD",
            "eqv2-omat": "EquiformerV2-31M-S2EF-OC20-All+MD",
            "eqv2-mptrj": "EquiformerV2-31M-S2EF-OC20-All+MD",
        }
        
        checkpoint = model_paths.get(variant, model_paths["base"])
        
        return OCPCalculator(
            checkpoint_path=checkpoint,
            cpu=self.device == "cpu"
        )
    
    def _load_sevennet(self, model_name: str, variant: str) -> Any:
        """加载 SevenNet 模型"""
        from sevenn.sevennet_calculator import SevenNetCalculator
        
        model_names = {
            "0": "SevenNet-0_11July2024",
            "mf-ompa": "SevenNet-MF-OMat24",
            "l3i5": "SevenNet-l3i5",
        }
        
        model = model_names.get(variant, model_names["0"])
        
        return SevenNetCalculator(
            model=model,
            device=self.device
        )


# 创建应用实例
worker = FAIRChemSevenNetWorker()
app = worker.app


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8002)
