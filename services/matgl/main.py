#!/usr/bin/env python3
"""
MatGL Worker Service (M3GNet + CHGNet)

支持模型:
- m3gnet, m3gnet-direct
- chgnet, chgnet-mpp
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import List, Any
from shared.worker_base import BaseWorkerApp


class MatGLWorker(BaseWorkerApp):
    """MatGL 模型服务 (M3GNet + CHGNet)"""
    
    MODELS = {
        "m3gnet": {"variant": "M3GNet-MP-2021.2.8-DIRECT-PES"},
        "m3gnet-direct": {"variant": "M3GNet-MP-2021.2.8-DIRECT-PES"},
        "chgnet": {"variant": "CHGNet-MPtrj-2024.2.13-PES-11M"},
        "chgnet-mpp": {"variant": "CHGNet-MPtrj-2024.2.13-PES-11M"},
    }
    
    def __init__(self):
        super().__init__(service_name="MatGL", port=8003)
    
    def get_available_models(self) -> List[str]:
        try:
            import matgl
            return list(self.MODELS.keys())
        except ImportError:
            print("MatGL not available")
            return []
    
    def load_model(self, model_name: str) -> Any:
        import matgl
        from matgl.ext.ase import M3GNetCalculator, CHGNetCalculator
        
        config = self.MODELS.get(model_name)
        if not config:
            raise ValueError(f"Unknown model: {model_name}")
        
        variant = config["variant"]
        
        if "M3GNet" in variant:
            pot = matgl.load_model(variant)
            return M3GNetCalculator(pot)
        elif "CHGNet" in variant:
            pot = matgl.load_model(variant)
            return CHGNetCalculator(pot)
        else:
            raise ValueError(f"Unknown variant: {variant}")


# 创建应用实例
worker = MatGLWorker()
app = worker.app


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8003)
