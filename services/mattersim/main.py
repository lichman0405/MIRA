#!/usr/bin/env python3
"""
MatterSim Worker Service

支持模型:
- mattersim-5m
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import List, Any
from shared.worker_base import BaseWorkerApp


class MatterSimWorker(BaseWorkerApp):
    """MatterSim 模型服务 (Microsoft)"""
    
    MODELS = {
        "mattersim-5m": {"variant": "MatterSim-v1.0.0-5M"},
    }
    
    def __init__(self):
        super().__init__(service_name="MatterSim", port=8005)
    
    def get_available_models(self) -> List[str]:
        try:
            from mattersim.forcefield import MatterSimCalculator
            return list(self.MODELS.keys())
        except ImportError:
            print("MatterSim not available")
            return []
    
    def load_model(self, model_name: str) -> Any:
        from mattersim.forcefield import MatterSimCalculator
        
        config = self.MODELS.get(model_name)
        if not config:
            raise ValueError(f"Unknown model: {model_name}")
        
        return MatterSimCalculator(device=self.device)


# 创建应用实例
worker = MatterSimWorker()
app = worker.app


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8005)
