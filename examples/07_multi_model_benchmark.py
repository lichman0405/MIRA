"""
MIRA Examples - 多模型基准测试
对所有可用模型进行系统性基准测试

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
import time
from pathlib import Path
from datetime import datetime
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


def run_benchmark_suite(
    atoms,
    model_names: List[str],
    run_optimization: bool = True,
    run_single_point: bool = True
) -> Dict[str, Any]:
    """
    运行完整的基准测试套件
    
    Args:
        atoms: ASE Atoms 对象
        model_names: 要测试的模型列表
        run_optimization: 是否运行优化测试
        run_single_point: 是否运行单点能量测试
        
    Returns:
        基准测试结果
    """
    atoms_data = atoms_to_dict(atoms)
    results = {
        "structure": atoms.get_chemical_formula(),
        "num_atoms": len(atoms),
        "timestamp": datetime.now().isoformat(),
        "models": {}
    }
    
    for model in model_names:
        print(f"\n{'='*60}")
        print(f"测试模型: {model}")
        print(f"{'='*60}")
        
        model_results = {}
        
        # 单点能量计算
        if run_single_point:
            print("\n[1] 单点能量计算...")
            start_time = time.time()
            
            response = requests.post(
                f"{BASE_URL}/single_point",
                json={
                    "atoms": atoms_data,
                    "model_name": model,
                    "compute_stress": True,
                    "compute_forces": True
                },
                timeout=300
            )
            
            elapsed = time.time() - start_time
            
            if response.ok:
                sp_result = response.json()
                model_results["single_point"] = {
                    "energy": sp_result.get("energy"),
                    "energy_per_atom": sp_result.get("energy_per_atom"),
                    "max_force": sp_result.get("max_force"),
                    "time_seconds": elapsed
                }
                print(f"    能量: {sp_result.get('energy', 0):.4f} eV")
                print(f"    耗时: {elapsed:.2f} s")
            else:
                model_results["single_point"] = {"error": response.text}
                print(f"    失败: {response.text[:100]}")
        
        # 结构优化
        if run_optimization:
            print("\n[2] 结构优化 (快速)...")
            start_time = time.time()
            
            response = requests.post(
                f"{BASE_URL}/optimization",
                json={
                    "atoms": atoms_data,
                    "model_name": model,
                    "fmax": 0.1,  # 宽松标准
                    "max_steps": 20,  # 少量步数
                    "optimizer": "BFGS",
                    "use_d3": True,
                    "fix_cell": True  # 固定晶胞加快速度
                },
                timeout=300
            )
            
            elapsed = time.time() - start_time
            
            if response.ok:
                opt_result = response.json()
                model_results["optimization"] = {
                    "initial_energy": opt_result.get("initial_energy"),
                    "final_energy": opt_result.get("final_energy"),
                    "steps": opt_result.get("steps"),
                    "converged": opt_result.get("converged"),
                    "time_seconds": elapsed
                }
                print(f"    初始能量: {opt_result.get('initial_energy', 0):.4f} eV")
                print(f"    最终能量: {opt_result.get('final_energy', 0):.4f} eV")
                print(f"    步数: {opt_result.get('steps', 0)}")
                print(f"    耗时: {elapsed:.2f} s")
            else:
                model_results["optimization"] = {"error": response.text}
                print(f"    失败: {response.text[:100]}")
        
        results["models"][model] = model_results
    
    return results


def print_benchmark_summary(results: Dict[str, Any]):
    """打印基准测试摘要"""
    print("\n" + "=" * 80)
    print("基准测试摘要")
    print("=" * 80)
    print(f"结构: {results['structure']}")
    print(f"原子数: {results['num_atoms']}")
    print(f"时间: {results['timestamp']}")
    print()
    
    # 单点能量比较
    print("单点能量比较:")
    print(f"{'模型':<25} {'能量 (eV)':>12} {'每原子 (eV)':>12} {'耗时 (s)':>10}")
    print("-" * 60)
    
    for model, data in results["models"].items():
        sp = data.get("single_point", {})
        if "error" not in sp:
            print(f"{model:<25} {sp.get('energy', 0):>12.4f} {sp.get('energy_per_atom', 0):>12.4f} {sp.get('time_seconds', 0):>10.2f}")
        else:
            print(f"{model:<25} {'错误':>12}")
    
    # 优化比较
    print("\n优化比较:")
    print(f"{'模型':<25} {'初始 (eV)':>12} {'最终 (eV)':>12} {'步数':>8} {'耗时 (s)':>10}")
    print("-" * 70)
    
    for model, data in results["models"].items():
        opt = data.get("optimization", {})
        if "error" not in opt:
            print(f"{model:<25} {opt.get('initial_energy', 0):>12.4f} {opt.get('final_energy', 0):>12.4f} {opt.get('steps', 0):>8} {opt.get('time_seconds', 0):>10.2f}")
        else:
            print(f"{model:<25} {'错误':>12}")


# ========== 主程序 ==========
if __name__ == "__main__":
    print("=" * 60)
    print("MIRA 多模型基准测试")
    print("=" * 60)
    
    # 获取可用模型
    available_models = get_service_models()
    if not available_models:
        print("[错误] 没有可用的模型")
        exit(1)
    
    print(f"\n可用模型 ({len(available_models)} 个):")
    for m in available_models:
        print(f"  - {m}")
    
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
    
    # 运行基准测试
    results = run_benchmark_suite(
        atoms,
        available_models,
        run_optimization=True,
        run_single_point=True
    )
    
    # 打印摘要
    print_benchmark_summary(results)
    
    # 保存结果
    output_file = Path(__file__).parent / "benchmark_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n结果已保存至: {output_file}")
    
    print("\n" + "=" * 60)
    print("基准测试完成！")
    print("=" * 60)
