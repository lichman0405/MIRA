"""
MIRA 配置管理
"""
from pydantic_settings import BaseSettings
from pydantic import Field
from typing import List, Optional
from pathlib import Path
from functools import lru_cache


class Settings(BaseSettings):
    """应用配置"""
    
    # 服务配置
    APP_NAME: str = "MIRA"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False
    HOST: str = "0.0.0.0"
    PORT: int = 8000
    
    # 路径配置
    BASE_DIR: Path = Path(__file__).parent.parent
    DATA_DIR: Path = Field(default_factory=lambda: Path("./data"))
    MODELS_CACHE_DIR: Path = Field(default_factory=lambda: Path("./models_cache"))
    RESULTS_DIR: Path = Field(default_factory=lambda: Path("./results"))
    STRUCTURES_DIR: Path = Field(default_factory=lambda: Path("./data/structures"))
    
    # 模型配置
    DEFAULT_DEVICE: str = "cuda"
    DEFAULT_DTYPE: str = "float32"
    PRELOAD_MODELS: List[str] = Field(default_factory=list)
    
    # 任务配置
    MAX_CONCURRENT_TASKS: int = 4
    TASK_TIMEOUT_SECONDS: int = 86400  # 24小时
    
    # 数据库配置（可选，用于持久化任务状态）
    DATABASE_URL: Optional[str] = None
    
    # Redis 配置（可选，用于任务队列）
    REDIS_URL: Optional[str] = None
    
    # CORS 配置
    CORS_ORIGINS: List[str] = Field(default_factory=lambda: ["*"])
    
    model_config = {
        "env_file": ".env",
        "env_file_encoding": "utf-8",
        "extra": "ignore"
    }
    
    def ensure_directories(self) -> None:
        """确保所有必要目录存在"""
        for dir_path in [self.DATA_DIR, self.MODELS_CACHE_DIR, 
                         self.RESULTS_DIR, self.STRUCTURES_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)


@lru_cache()
def get_settings() -> Settings:
    """获取配置单例"""
    settings = Settings()
    settings.ensure_directories()
    return settings


# 导出配置实例
settings = get_settings()
