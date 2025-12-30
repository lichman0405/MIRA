"""
MIRA 客户端工具

提供统一的服务检查和配置功能。
使用 Docker 微服务时无需本地安装 ML 模型包。
"""
import os
import requests
from typing import List, Optional


def get_gateway_url() -> str:
    """获取网关 URL"""
    return os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000")


def get_api_base_url() -> str:
    """获取 API 基础 URL"""
    return f"{get_gateway_url()}/api/v1"


def check_service_available(timeout: float = 3.0) -> bool:
    """
    检查 MIRA 服务是否可用
    
    Args:
        timeout: 请求超时时间（秒）
    
    Returns:
        服务是否可用
    """
    try:
        # 健康检查端点在根路径，不在 /api/v1 下
        resp = requests.get(f"{get_gateway_url()}/health", timeout=timeout)
        return resp.ok
    except Exception:
        return False


def get_service_models() -> List[str]:
    """
    从服务获取可用模型列表
    
    Returns:
        可用模型的 name 列表
    """
    try:
        resp = requests.get(f"{get_api_base_url()}/models/", timeout=10)
        if resp.ok:
            data = resp.json()
            # Gateway API 返回格式: {"total": N, "models": [...]}
            models = data.get("models", []) if isinstance(data, dict) else data
            return [m["name"] for m in models]
    except Exception:
        pass
    return []


def ensure_service_ready(verbose: bool = True) -> bool:
    """
    确保服务可用，提供用户友好的提示
    
    Args:
        verbose: 是否打印详细信息
    
    Returns:
        服务是否可用
    """
    gateway_url = get_gateway_url()
    
    if check_service_available():
        if verbose:
            models = get_service_models()
            print(f"[✓] 服务已连接: {gateway_url}")
            print(f"[✓] 可用模型: {len(models)} 个")
        return True
    
    # 服务不可用
    if verbose:
        print(f"[✗] 服务不可用: {gateway_url}")
        print()
        print("请先启动 MIRA 服务:")
        print("  Docker (推荐): ./scripts/deploy.sh test-cpu")
        print("  GPU Docker:    ./scripts/deploy.sh test")
        print()
        print("或设置远程服务器地址:")
        print("  export MIRA_GATEWAY_URL=http://192.168.100.207:8000")
    
    return False


def init_client(verbose: bool = True, exit_on_fail: bool = True):
    """
    初始化客户端，检查服务并返回配置
    
    Args:
        verbose: 是否打印详细信息
        exit_on_fail: 服务不可用时是否退出
    
    Returns:
        tuple: (gateway_url, base_url)
    """
    gateway_url = get_gateway_url()
    base_url = get_api_base_url()
    
    if not ensure_service_ready(verbose=verbose):
        if exit_on_fail:
            import sys
            sys.exit(1)
    
    return gateway_url, base_url
