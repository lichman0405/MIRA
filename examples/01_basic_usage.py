"""
MIRA Examples - 基础使用示例
演示如何连接服务、查询模型、查看服务状态

运行前确保:
1. 服务已启动: 
   - Docker (推荐): ./scripts/deploy.sh test-cpu
   - GPU Docker: ./scripts/deploy.sh test

支持环境变量:
    MIRA_GATEWAY_URL=http://192.168.100.207:8000  # 测试服务器

注意: 使用 Docker 微服务时，无需本地安装 ML 模型包
"""
import requests
from pathlib import Path

# 客户端初始化
from client_utils import init_client

GATEWAY_URL, BASE_URL = init_client(verbose=True)


# ========== 1. 健康检查 ==========
def check_health():
    """检查服务是否正常运行，显示 Worker 状态"""
    response = requests.get(f"{GATEWAY_URL}/health")
    data = response.json()
    
    print("\n=== 服务健康检查 ===")
    print(f"Gateway 状态: {data.get('gateway', 'unknown')}")
    print(f"可用 Workers: {data.get('available_workers', 0)}/{data.get('total_workers', 0)}")
    
    workers = data.get("workers", {})
    print("\nWorker 详情:")
    for name, info in workers.items():
        status = info.get("status", "unknown")
        if status == "healthy":
            device = info.get("device", "unknown")
            models = info.get("models", [])
            print(f"  ✓ {name}: {len(models)} 个模型 (设备: {device})")
        else:
            error = info.get("error", "未知错误")
            print(f"  ✗ {name}: {error[:50]}...")
    
    return data.get("status") == "healthy"


# ========== 2. 列出所有可用模型 ==========
def list_all_models():
    """获取所有可用的 ML 力场模型"""
    response = requests.get(f"{BASE_URL}/models/")
    data = response.json()
    
    # Gateway API 返回格式: {"total": N, "models": [...]}
    models = data.get("models", []) if isinstance(data, dict) else data
    
    print("\n=== 可用模型列表 ===")
    print(f"共 {len(models)} 个模型\n")
    
    # 按服务分组显示
    services = {}
    for m in models:
        # Gateway 格式: {"name": "mace-mp", "service": "mace-orb", "url": "..."}
        service = m.get('service', 'unknown')
        if service not in services:
            services[service] = []
        services[service].append(m)
    
    for service, model_list in services.items():
        print(f"【{service}】")
        for m in model_list:
            print(f"  - {m['name']}")
        print()
    
    return models


# ========== 3. 列出服务状态 ==========
def list_services():
    """获取所有服务的详细状态"""
    response = requests.get(f"{BASE_URL}/services")
    data = response.json()
    
    print("\n=== 服务状态 ===")
    for name, info in data.items():
        status = info.get("status", "unknown")
        if status == "healthy":
            gpu = info.get("gpu_name", "N/A")
            device = info.get("device", "unknown")
            print(f"  ✓ {name}")
            print(f"      服务: {info.get('service', 'N/A')}")
            print(f"      设备: {device} ({gpu})")
            print(f"      模型: {', '.join(info.get('models', []))}")
        else:
            print(f"  ✗ {name}: 不可用")
    
    return data


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 基础使用示例")
    print("=" * 60)
    
    # 1. 健康检查
    if not check_health():
        print("\n[警告] 部分服务不可用，但可以继续使用已启动的 Worker")
    
    # 2. 列出可用模型
    models = list_all_models()
    
    if not models:
        print("\n[警告] 没有可用的模型，请检查 Worker 服务")
        exit(1)
    
    # 3. 列出服务状态
    list_services()
    
    print("\n" + "=" * 60)
    print("基础示例完成！")
    print("=" * 60)
    print("\n下一步:")
    print("  - 运行 02_structure_optimization.py 进行结构优化")
    print("  - 运行 03_md_stability.py 进行 MD 稳定性测试")
