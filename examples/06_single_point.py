"""
MIRA Examples - 单点能量计算
计算结构的单点能量、力和应力

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
from typing import Dict, Any, List
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


def run_single_point(
    atoms,
    model_name: str,
    compute_stress: bool = True,
    compute_forces: bool = True
) -> Dict[str, Any]:
    """
    计算单点能量
    
    Args:
        atoms: ASE Atoms 对象
        model_name: 模型名称
        compute_stress: 是否计算应力
        compute_forces: 是否计算力
        
    Returns:
        计算结果
    """
    if hasattr(atoms, 'get_chemical_symbols'):
        atoms_data = atoms_to_dict(atoms)
    else:
        atoms_data = atoms
    
    request_data = {
        "atoms": atoms_data,
        "model_name": model_name,
        "compute_stress": compute_stress,
        "compute_forces": compute_forces
    }
    
    print(f"\n=== 单点能量计算 ===")
    print(f"  模型: {model_name}")
    print("\n正在计算...")
    
    response = requests.post(
        f"{BASE_URL}/single_point",
        json=request_data,
        timeout=300
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 计算完成 ===")
        print(f"  总能量: {result.get('energy', 'N/A'):.6f} eV")
        print(f"  每原子能量: {result.get('energy_per_atom', 'N/A'):.6f} eV/atom")
        
        if 'max_force' in result:
            print(f"  最大力: {result['max_force']:.6f} eV/Å")
        
        if 'stress' in result:
            stress = result['stress']
            print(f"  应力张量 (GPa): {stress}")
        
        return result
    else:
        print(f"\n计算失败: {response.text}")
        return {"error": response.text}


def compare_models_single_point(
    atoms,
    model_names: List[str]
) -> List[Dict[str, Any]]:
    """
    使用多个模型计算单点能量进行比较
    """
    results = []
    
    print(f"\n{'='*60}")
    print(f"多模型单点能量比较")
    print(f"结构: {atoms.get_chemical_formula()}")
    print(f"原子数: {len(atoms)}")
    print(f"{'='*60}")
    
    for model in model_names:
        print(f"\n>>> 模型: {model}")
        result = run_single_point(atoms, model)
        result["model"] = model
        results.append(result)
    
    # 总结
    print(f"\n{'='*60}")
    print("计算结果汇总")
    print(f"{'='*60}")
    print(f"{'模型':<20} {'总能量 (eV)':>15} {'每原子 (eV)':>12} {'最大力 (eV/Å)':>15}")
    print("-" * 65)
    
    for r in results:
        if "error" not in r:
            print(f"{r['model']:<20} {r.get('energy', 0):>15.4f} {r.get('energy_per_atom', 0):>12.4f} {r.get('max_force', 0):>15.4f}")
        else:
            print(f"{r['model']:<20} {'错误':>15}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 单点能量计算示例")
    print("=" * 60)
    
    # 获取可用模型
    available_models = get_service_models()
    if not available_models:
        print("[错误] 没有可用的模型")
        exit(1)
    
    print(f"\n可用模型: {', '.join(available_models)}")
    
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
    
    # 单模型计算（避免连续请求导致 Worker 崩溃）
    model_to_use = available_models[0]
    print(f"\n使用模型: {model_to_use}")
    run_single_point(atoms, model_to_use)
    
    print("\n" + "=" * 60)
    print("单点能量计算完成！")
    print("=" * 60)
