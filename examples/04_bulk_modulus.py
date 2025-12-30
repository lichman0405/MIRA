"""
MIRA Examples - 体积模量计算
通过 E-V 曲线拟合 Birch-Murnaghan 状态方程计算体积模量

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
from typing import Dict, Any
from ase.io import read

# 客户端初始化
from client_utils import init_client, get_service_models


GATEWAY_URL, BASE_URL = init_client(verbose=True)


def atoms_to_dict(atoms) -> Dict[str, Any]:
    """将 ASE Atoms 对象转换为 API 请求格式"""
    return {
        "symbols": list(atoms.get_chemical_symbols()),
        "positions": atoms.get_positions().tolist(),
        "cell": atoms.get_cell().tolist() if atoms.cell is not None else None,
        "pbc": [bool(p) for p in atoms.get_pbc()]  # 转换 numpy bool 为 Python bool
    }


def run_bulk_modulus(
    atoms,
    model_name: str,
    strain_range: float = 0.06,
    num_points: int = 7,
    use_d3: bool = True
) -> Dict[str, Any]:
    """
    计算体积模量
    
    Args:
        atoms: ASE Atoms 对象
        model_name: 模型名称
        strain_range: 应变范围 (±)
        num_points: 采样点数
        use_d3: 是否启用 D3 校正
        
    Returns:
        体积模量计算结果
    """
    if hasattr(atoms, 'get_chemical_symbols'):
        atoms_data = atoms_to_dict(atoms)
    else:
        atoms_data = atoms
    
    request_data = {
        "atoms": atoms_data,
        "model_name": model_name,
        "strain_range": strain_range,
        "num_points": num_points,
        "use_d3": use_d3
    }
    
    print(f"\n=== 体积模量计算 ===")
    print(f"  模型: {model_name}")
    print(f"  应变范围: ±{strain_range*100:.1f}%")
    print(f"  采样点数: {num_points}")
    print("\n正在计算 E-V 曲线...")
    
    response = requests.post(
        f"{BASE_URL}/bulk_modulus",
        json=request_data,
        timeout=600
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 计算完成 ===")
        print(f"  体积模量 B0: {result.get('bulk_modulus', 'N/A'):.2f} GPa")
        print(f"  平衡体积 V0: {result.get('equilibrium_volume', 'N/A'):.2f} Å³")
        print(f"  B0 导数 B': {result.get('bulk_modulus_derivative', 'N/A'):.2f}")
        print(f"  拟合 R²: {result.get('r_squared', 'N/A'):.4f}")
        return result
    else:
        print(f"\n计算失败: {response.text}")
        return {"error": response.text}


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 体积模量计算示例")
    print("=" * 60)
    
    # 获取可用模型
    available_models = get_service_models()
    if not available_models:
        print("[错误] 没有可用的模型")
        exit(1)
    
    print(f"\n可用模型: {', '.join(available_models)}")
    model_to_use = available_models[0]
    print(f"将使用模型: {model_to_use}")
    
    # 加载示例结构 - 优先选择小结构
    structures_dir = Path(__file__).parent / "structures"
    preferred_files = ["ZIF-8.cif", "HKUST-1.cif", "UiO-66.cif", "MOF-5.cif"]
    structure_file = None
    
    for name in preferred_files:
        candidate = structures_dir / name
        if candidate.exists():
            structure_file = candidate
            break
    
    if structure_file is None:
        cif_files = list(structures_dir.glob("*.cif"))
        if cif_files:
            structure_file = cif_files[0]
    
    if structure_file is None:
        print(f"\n[警告] 在 {structures_dir} 中没有找到 CIF 文件")
        exit(1)
    
    print(f"\n加载结构: {structure_file.name}")
    
    try:
        atoms = read(str(structure_file))
        print(f"  化学式: {atoms.get_chemical_formula()}")
        print(f"  原子数: {len(atoms)}")
    except Exception as e:
        print(f"加载结构失败: {e}")
        exit(1)
    
    # 计算体积模量
    result = run_bulk_modulus(
        atoms,
        model_name=model_to_use,
        strain_range=0.04,  # ±4%
        num_points=5,
        use_d3=True
    )
    
    print("\n" + "=" * 60)
    print("体积模量计算完成！")
    print("=" * 60)
