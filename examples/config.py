"""
MIRA 客户端配置

支持通过环境变量或配置文件设置服务地址
"""
import os
from dataclasses import dataclass
from typing import Optional


@dataclass
class MIRAConfig:
    """MIRA 服务配置"""
    
    # 服务地址
    gateway_url: str = "http://localhost:8000"
    
    # 超时设置 (秒)
    timeout: float = 600.0
    
    # 是否显示详细信息
    verbose: bool = True
    
    @classmethod
    def from_env(cls) -> "MIRAConfig":
        """从环境变量加载配置"""
        return cls(
            gateway_url=os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000"),
            timeout=float(os.getenv("MIRA_TIMEOUT", "600")),
            verbose=os.getenv("MIRA_VERBOSE", "true").lower() == "true"
        )
    
    @classmethod
    def for_server(cls, ip: str, port: int = 8000) -> "MIRAConfig":
        """为特定服务器创建配置"""
        return cls(gateway_url=f"http://{ip}:{port}")
    
    @property
    def api_base_url(self) -> str:
        """API 基础 URL"""
        return f"{self.gateway_url}/api/v1"


# ============================================
# 预定义配置
# ============================================

# 本地开发环境
LOCAL_CONFIG = MIRAConfig(gateway_url="http://localhost:8000")

# 测试服务器 (192.168.100.207)
TEST_SERVER_CONFIG = MIRAConfig(gateway_url="http://192.168.100.207:8000")

# 从环境变量加载
ENV_CONFIG = MIRAConfig.from_env()


def get_config(server: Optional[str] = None) -> MIRAConfig:
    """
    获取配置
    
    Args:
        server: 服务器标识 ("local", "test", None=从环境变量)
    
    Returns:
        MIRAConfig 配置对象
    """
    if server == "local":
        return LOCAL_CONFIG
    elif server == "test":
        return TEST_SERVER_CONFIG
    else:
        return ENV_CONFIG


# 便捷函数
def get_base_url(server: Optional[str] = None) -> str:
    """获取 API 基础 URL"""
    return get_config(server).api_base_url
