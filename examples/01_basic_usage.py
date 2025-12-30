"""
MIRA Examples - 基础使用示例
演示如何连接服务、查询模型、上传结构

运行前确保:
1. 已安装依赖: python examples/setup_check.py
2. 服务已启动: uvicorn app.main:app --host 0.0.0.0 --port 8000
"""
import requests
from pathlib import Path

# 依赖检查
try:
    from setup_check import ensure_dependencies, check_available_models
    ensure_dependencies(verbose=False)
except ImportError:
    print("提示: 运行 'python examples/setup_check.py' 检查依赖")

# ========== 配置 ==========
BASE_URL = "http://localhost:8000/api/v1"

# ========== 1. 健康检查 ==========
def check_health():
    """检查服务是否正常运行"""
    response = requests.get(f"{BASE_URL}/health")
    print("=== 健康检查 ===")
    print(response.json())
    return response.ok


# ========== 2. 列出所有可用模型 ==========
def list_all_models():
    """获取所有可用的 ML 力场模型"""
    response = requests.get(f"{BASE_URL}/models/")
    models = response.json()
    
    print("\n=== 可用模型列表 ===")
    print(f"共 {len(models)} 个模型\n")
    
    # 按家族分组显示
    families = {}
    for m in models:
        family = m['family']
        if family not in families:
            families[family] = []
        families[family].append(m)
    
    for family, model_list in families.items():
        print(f"【{family}】")
        for m in model_list:
            print(f"  - {m['key']}: {m['name']}")
        print()
    
    return models


# ========== 3. 列出模型家族 ==========
def list_model_families():
    """获取所有模型家族"""
    response = requests.get(f"{BASE_URL}/models/families")
    families = response.json()
    
    print("\n=== 模型家族 ===")
    for f in families:
        print(f"  - {f}")
    
    return families


# ========== 4. 获取特定模型信息 ==========
def get_model_info(model_key: str):
    """获取特定模型的详细信息"""
    response = requests.get(f"{BASE_URL}/models/{model_key}")
    
    if response.ok:
        info = response.json()
        print(f"\n=== 模型信息: {model_key} ===")
        print(f"名称: {info['name']}")
        print(f"家族: {info['family']}")
        print(f"描述: {info.get('description', 'N/A')}")
        return info
    else:
        print(f"模型 {model_key} 未找到")
        return None


# ========== 5. 上传结构文件 ==========
def upload_structure(file_path: str, name: str = None):
    """上传 CIF/XYZ 结构文件"""
    path = Path(file_path)
    
    if not path.exists():
        print(f"文件不存在: {file_path}")
        return None
    
    # 确定格式
    ext = path.suffix.lower()
    format_map = {'.cif': 'cif', '.xyz': 'xyz', '.vasp': 'poscar'}
    file_format = format_map.get(ext, 'cif')
    
    # 读取内容
    content = path.read_text()
    
    # 上传
    response = requests.post(
        f"{BASE_URL}/structures/upload",
        data={
            "name": name or path.stem,
            "format": file_format,
            "content": content
        }
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 结构上传成功 ===")
        print(f"ID: {result['id']}")
        print(f"名称: {result['name']}")
        print(f"化学式: {result['formula']}")
        print(f"原子数: {result['num_atoms']}")
        print(f"体积: {result['cell_volume']:.2f} Å³")
        return result
    else:
        print(f"上传失败: {response.text}")
        return None


# ========== 6. 列出已上传的结构 ==========
def list_structures():
    """列出所有已上传的结构"""
    response = requests.get(f"{BASE_URL}/structures/")
    
    if response.ok:
        data = response.json()
        print(f"\n=== 已上传结构 ({data['total']} 个) ===")
        for s in data['structures']:
            print(f"  [{s['id'][:8]}...] {s['name']} - {s['formula']} ({s['num_atoms']} atoms)")
        return data['structures']
    
    return []


# ========== 7. 获取结构预览数据 ==========
def get_structure_preview(structure_id: str):
    """获取结构的坐标数据（用于可视化）"""
    response = requests.get(f"{BASE_URL}/structures/{structure_id}/preview")
    
    if response.ok:
        preview = response.json()
        print(f"\n=== 结构预览: {preview['name']} ===")
        print(f"化学式: {preview['formula']}")
        print(f"原子数: {len(preview['symbols'])}")
        print(f"元素: {set(preview['symbols'])}")
        print(f"PBC: {preview['pbc']}")
        return preview
    
    return None


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 基础使用示例")
    print("=" * 60)
    
    # 1. 健康检查
    if not check_health():
        print("服务未运行，请先启动 MIRA 服务")
        exit(1)
    
    # 2. 列出模型
    list_all_models()
    
    # 3. 获取特定模型信息
    get_model_info("mace-mp")
    get_model_info("orb-v2")
    
    # 4. 上传示例结构
    structures_dir = Path(__file__).parent / "structures"
    
    for cif_file in structures_dir.glob("*.cif"):
        upload_structure(str(cif_file))
    
    # 5. 列出已上传结构
    structures = list_structures()
    
    # 6. 获取第一个结构的预览
    if structures:
        get_structure_preview(structures[0]['id'])
    
    print("\n" + "=" * 60)
    print("基础示例完成！")
    print("=" * 60)
