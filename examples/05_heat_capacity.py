"""
MIRA Examples - 热容计算
使用 Phonopy 进行声子计算，获得热容和其他热力学性质

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


def run_heat_capacity(
    atoms,
    model_name: str,
    temperatures: List[float] = None,
    supercell: List[int] = None,
    use_d3: bool = True
) -> Dict[str, Any]:
    """
    计算热容
    
    Args:
        atoms: ASE Atoms 对象
        model_name: 模型名称
        temperatures: 温度列表 (K)
        supercell: 超胞尺寸 [nx, ny, nz]
        use_d3: 是否启用 D3 校正
        
    Returns:
        热容计算结果
    """
    if temperatures is None:
        temperatures = [100, 200, 300, 400, 500]
    if supercell is None:
        supercell = [2, 2, 2]
    
    if hasattr(atoms, 'get_chemical_symbols'):
        atoms_data = atoms_to_dict(atoms)
    else:
        atoms_data = atoms
    
    request_data = {
        "atoms": atoms_data,
        "model_name": model_name,
        "temperatures": temperatures,
        "supercell": supercell,
        "use_d3": use_d3
    }
    
    print(f"\n=== 热容计算 ===")
    print(f"  模型: {model_name}")
    print(f"  温度范围: {min(temperatures)}-{max(temperatures)} K")
    print(f"  超胞尺寸: {supercell}")
    print("\n正在计算声子...")
    
    response = requests.post(
        f"{BASE_URL}/heat_capacity",
        json=request_data,
        timeout=1800  # 30 分钟超时
    )
    
    if response.ok:
        result = response.json()
        print(f"\n=== 计算完成 ===")
        
        # 显示各温度下的热容
        cv_values = result.get('heat_capacity', [])
        temps = result.get('temperatures', temperatures)
        
        print(f"\n温度 (K)  |  Cv (J/mol·K)")
        print("-" * 30)
        for t, cv in zip(temps, cv_values):
            print(f"  {t:>6.1f}  |  {cv:>8.2f}")
        
        if 'has_imaginary' in result:
            print(f"\n虚频模式: {'存在' if result['has_imaginary'] else '无'}")
        
        return result
    else:
        print(f"\n计算失败: {response.text}")
        return {"error": response.text}


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 热容计算示例")
    print("=" * 60)
    
    # 获取可用模型
    available_models = get_service_models()
    if not available_models:
        print("[错误] 没有可用的模型")
        exit(1)
    
    print(f"\n可用模型: {', '.join(available_models)}")
    model_to_use = available_models[0]
    print(f"将使用模型: {model_to_use}")
    
    # 加载示例结构
    structures_dir = Path(__file__).parent / "structures"
    cif_files = list(structures_dir.glob("*.cif"))
    
    if not cif_files:
        print(f"\n[警告] 在 {structures_dir} 中没有找到 CIF 文件")
        exit(1)
    
    structure_file = cif_files[0]
    print(f"\n加载结构: {structure_file.name}")
    
    try:
        atoms = read(str(structure_file))
        print(f"  化学式: {atoms.get_chemical_formula()}")
        print(f"  原子数: {len(atoms)}")
    except Exception as e:
        print(f"加载结构失败: {e}")
        exit(1)
    
    # 热容计算（声子计算比较慢，使用小的超胞）
    print("\n[提示] 声子计算可能需要较长时间...")
    
    result = run_heat_capacity(
        atoms,
        model_name=model_to_use,
        temperatures=[100, 200, 300],
        supercell=[1, 1, 1],  # 小超胞用于快速测试
        use_d3=True
    )
    
    print("\n" + "=" * 60)
    print("热容计算完成！")
    print("=" * 60)
