"""
MIRA Examples - MD 稳定性测试
对 MOF 结构进行分子动力学模拟，评估热稳定性

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


def run_stability_test(
    atoms,
    model_name: str,
    temperature: float = 300.0,
    pressure: float = 0.0,
    timestep: float = 1.0,
    equilibration_steps: int = 500,
    production_steps: int = 1000,
    use_d3: bool = True
) -> Dict[str, Any]:
    """
    运行 MD 稳定性测试
    
    Args:
        atoms: ASE Atoms 对象
        model_name: 模型名称
        temperature: 温度 (K)
        pressure: 压力 (GPa), 0 表示 NVT
        timestep: 时间步长 (fs)
        equilibration_steps: 平衡步数
        production_steps: 生产步数
        use_d3: 是否启用 D3 校正
        
    Returns:
        稳定性测试结果
    """
    if hasattr(atoms, 'get_chemical_symbols'):
        atoms_data = atoms_to_dict(atoms)
    else:
        atoms_data = atoms
    
    request_data = {
        "atoms": atoms_data,
        "model_name": model_name,
        "temperature": temperature,
        "pressure": pressure,
        "timestep": timestep,
        "equilibration_steps": equilibration_steps,
        "production_steps": production_steps,
        "use_d3": use_d3
    }
    
    print(f"\n=== MD 稳定性测试 ===")
    print(f"  模型: {model_name}")
    print(f"  温度: {temperature} K")
    print(f"  压力: {pressure} GPa")
    print(f"  时间步长: {timestep} fs")
    print(f"  平衡步数: {equilibration_steps}")
    print(f"  生产步数: {production_steps}")
    print("\n正在运行 MD 模拟...")
    
    response = requests.post(
        f"{BASE_URL}/stability",
        json=request_data,
        timeout=1800  # 30 分钟超时
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 测试完成 ===")
        print(f"  最终温度: {result.get('final_temperature', 'N/A'):.1f} K")
        print(f"  能量漂移: {result.get('energy_drift', 'N/A'):.6f} eV/atom")
        print(f"  结构稳定: {'是' if result.get('is_stable', False) else '否'}")
        if 'rmsd' in result:
            print(f"  RMSD: {result['rmsd']:.4f} Å")
        return result
    else:
        print(f"\n测试失败: {response.text}")
        return {"error": response.text}


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA MD 稳定性测试示例")
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
    
    # 运行稳定性测试（使用较少步数进行快速测试）
    result = run_stability_test(
        atoms,
        model_name=model_to_use,
        temperature=300.0,
        equilibration_steps=100,  # 快速测试
        production_steps=200,
        use_d3=True
    )
    
    print("\n" + "=" * 60)
    print("MD 稳定性测试完成！")
    print("=" * 60)
