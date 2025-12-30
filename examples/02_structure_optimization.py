"""
MIRA Examples - 结构优化示例
使用 ML 力场模型对 MOF 结构进行优化

运行前确保:
1. 服务已启动: 
   - Docker (推荐): ./scripts/deploy.sh test-cpu
   - GPU Docker: ./scripts/deploy.sh test

支持环境变量:
    MIRA_GATEWAY_URL=http://192.168.100.207:8000  # 测试服务器

注意: 使用 Docker 微服务时，无需本地安装 ML 模型包
"""
import requests
import json
from pathlib import Path
from typing import Dict, Any, List, Optional
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


def load_structure(file_path: str):
    """从文件加载结构"""
    return read(file_path)


def run_optimization(
    atoms,
    model_name: str,
    fmax: float = 0.05,
    max_steps: int = 200,
    optimizer: str = "BFGS",
    use_d3: bool = True,
    fix_cell: bool = False
) -> Dict[str, Any]:
    """
    对结构运行优化
    
    Args:
        atoms: ASE Atoms 对象或原子数据字典
        model_name: 模型名称 (如 "mace-mp", "orb-v2")
        fmax: 力收敛阈值 (eV/Å)
        max_steps: 最大优化步数
        optimizer: 优化器 (BFGS/FIRE/LBFGS)
        use_d3: 是否启用 D3 色散校正
        fix_cell: 是否固定晶胞
        
    Returns:
        优化结果字典
    """
    # 准备请求数据
    if hasattr(atoms, 'get_chemical_symbols'):
        atoms_data = atoms_to_dict(atoms)
    else:
        atoms_data = atoms
    
    request_data = {
        "atoms": atoms_data,
        "model_name": model_name,
        "fmax": fmax,
        "max_steps": max_steps,
        "optimizer": optimizer,
        "use_d3": use_d3,
        "fix_cell": fix_cell
    }
    
    print(f"\n优化参数:")
    print(f"  模型: {model_name}")
    print(f"  fmax: {fmax} eV/Å")
    print(f"  最大步数: {max_steps}")
    print(f"  优化器: {optimizer}")
    print(f"  D3 校正: {'启用' if use_d3 else '禁用'}")
    print(f"  固定晶胞: {'是' if fix_cell else '否'}")
    print("\n正在优化...")
    
    response = requests.post(
        f"{BASE_URL}/optimization",
        json=request_data,
        timeout=600  # 10 分钟超时
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 优化完成 ===")
        print(f"  初始能量: {result.get('initial_energy', 'N/A'):.6f} eV")
        print(f"  最终能量: {result.get('final_energy', 'N/A'):.6f} eV")
        print(f"  能量变化: {result.get('energy_change', 'N/A'):.6f} eV")
        print(f"  优化步数: {result.get('steps', 'N/A')}")
        print(f"  是否收敛: {'是' if result.get('converged', False) else '否'}")
        return result
    else:
        print(f"\n优化失败: {response.text}")
        return {"error": response.text}


def compare_models(
    atoms,
    model_names: List[str],
    fmax: float = 0.05,
    max_steps: int = 100
) -> List[Dict[str, Any]]:
    """
    使用多个模型对同一结构进行优化比较
    
    Args:
        atoms: ASE Atoms 对象
        model_names: 模型名称列表
        fmax: 力收敛阈值
        max_steps: 最大优化步数
        
    Returns:
        各模型的优化结果列表
    """
    results = []
    
    print(f"\n{'='*60}")
    print(f"多模型优化比较")
    print(f"结构: {atoms.get_chemical_formula()}")
    print(f"原子数: {len(atoms)}")
    print(f"{'='*60}")
    
    for model in model_names:
        print(f"\n>>> 模型: {model}")
        result = run_optimization(
            atoms, 
            model,
            fmax=fmax,
            max_steps=max_steps
        )
        result["model"] = model
        results.append(result)
    
    # 总结
    print(f"\n{'='*60}")
    print("优化结果汇总")
    print(f"{'='*60}")
    print(f"{'模型':<20} {'初始能量':>12} {'最终能量':>12} {'步数':>6} {'收敛':>6}")
    print("-" * 60)
    
    for r in results:
        if "error" not in r:
            converged = "✓" if r.get("converged") else "✗"
            print(f"{r['model']:<20} {r.get('initial_energy', 0):>12.4f} {r.get('final_energy', 0):>12.4f} {r.get('steps', 0):>6} {converged:>6}")
        else:
            print(f"{r['model']:<20} {'错误':>12}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 结构优化示例")
    print("=" * 60)
    
    # 获取可用模型
    available_models = get_service_models()
    if not available_models:
        print("[错误] 没有可用的模型")
        exit(1)
    
    print(f"\n可用模型: {', '.join(available_models)}")
    
    # 选择要使用的模型（使用第一个可用的）
    model_to_use = available_models[0]
    print(f"将使用模型: {model_to_use}")
    
    # 加载示例结构
    structures_dir = Path(__file__).parent / "structures"
    cif_files = list(structures_dir.glob("*.cif"))
    
    if not cif_files:
        print(f"\n[警告] 在 {structures_dir} 中没有找到 CIF 文件")
        print("请添加 MOF 结构文件后重试")
        exit(1)
    
    # 使用第一个结构文件
    structure_file = cif_files[0]
    print(f"\n加载结构: {structure_file.name}")
    
    try:
        atoms = load_structure(str(structure_file))
        print(f"  化学式: {atoms.get_chemical_formula()}")
        print(f"  原子数: {len(atoms)}")
    except Exception as e:
        print(f"加载结构失败: {e}")
        exit(1)
    
    # 运行单模型优化
    result = run_optimization(
        atoms,
        model_name=model_to_use,
        fmax=0.05,
        max_steps=100,
        use_d3=True
    )
    
    # 如果有多个可用模型，进行比较
    if len(available_models) >= 2:
        print("\n" + "=" * 60)
        print("多模型比较 (使用前2个模型)")
        print("=" * 60)
        
        compare_models(
            atoms,
            available_models[:2],
            fmax=0.1,  # 使用较宽松的收敛标准加快测试
            max_steps=50
        )
    
    print("\n" + "=" * 60)
    print("优化示例完成！")
    print("=" * 60)
