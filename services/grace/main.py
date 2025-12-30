#!/usr/bin/env python3
"""
GRACE Worker Service

支持模型:
- grace-2l, grace-2l-omat, grace-2m

注意: GRACE 使用 TensorFlow，与 PyTorch 模型隔离
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import List, Any
from shared.worker_base import BaseWorkerApp


class GRACEWorker(BaseWorkerApp):
    """GRACE 模型服务"""
    
    MODELS = {
        "grace-2l": {"variant": "grace-2l"},
        "grace-2l-omat": {"variant": "grace-2l-omat"},
        "grace-2m": {"variant": "grace-2m"},
    }
    
    def __init__(self):
        # GRACE 使用 CPU 或特定 GPU 配置
        super().__init__(service_name="GRACE", port=8004)
        # 限制 TensorFlow GPU 内存使用
        self._configure_tensorflow()
    
    def _configure_tensorflow(self):
        """配置 TensorFlow GPU 内存"""
        try:
            import tensorflow as tf
            gpus = tf.config.list_physical_devices('GPU')
            if gpus:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
        except Exception as e:
            print(f"TensorFlow GPU config failed: {e}")
    
    def get_available_models(self) -> List[str]:
        try:
            import tensorpotential
            return list(self.MODELS.keys())
        except ImportError:
            print("TensorPotential (GRACE) not available")
            return []
    
    def load_model(self, model_name: str) -> Any:
        from tensorpotential.calculator import TPCalculator
        
        config = self.MODELS.get(model_name)
        if not config:
            raise ValueError(f"Unknown model: {model_name}")
        
        variant = config["variant"]
        
        # GRACE 模型路径
        model_paths = {
            "grace-2l": "GRACE-2L-MP",
            "grace-2l-omat": "GRACE-2L-OMat24",
            "grace-2m": "GRACE-2M-MP",
        }
        
        model_name_full = model_paths.get(variant, model_paths["grace-2l"])
        
        return TPCalculator(model=model_name_full)


# 创建应用实例
worker = GRACEWorker()
app = worker.app


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8004)
