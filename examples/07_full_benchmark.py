"""
MIRA Examples - 全模型基准测试
对所有可用模型进行系统性基准测试

运行前确保:
1. 服务已启动: 
   - Docker (推荐): ./scripts/deploy.sh test-cpu
   - GPU Docker: ./scripts/deploy.sh test

支持环境变量:
    MIRA_GATEWAY_URL=http://192.168.100.207:8000  # 测试服务器

注意: 使用 Docker 微服务时，无需本地安装 ML 模型包
"""
import os
import requests
import time
import json
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List

# 客户端初始化
from client_utils import init_client, get_service_models

GATEWAY_URL, BASE_URL = init_client(verbose=True)

# 所有可用模型
ALL_MODELS = [
    # MACE 系列
    "mace-mp", "mace-off23", "mace-omat", "mace-mpa", "mace-ani",
    # ORB 系列
    "orb-v2", "orb-d3-v2", "orb-v3", "orb-v3-mpa", "orb-omat-v3-lora",
    # OMAT24 系列
    "omat24-base", "omat24-large", "eqv2-omat", "eqv2-mptrj",
    # GRACE 系列
    "grace-2l", "grace-2l-omat", "grace-2m",
    # 其他
    "mattersim-5m",
    "sevennet-0", "sevennet-mf-ompa", "sevennet-l3i5",
    "posegnn",
    "m3gnet", "chgnet"
]

# 快速测试子集
QUICK_TEST_MODELS = [
    "mace-mp",
    "orb-v2",
    "omat24-base",
    "grace-2l",
    "mattersim-5m",
    "sevennet-0"
]


def upload_structure(file_path: str) -> Optional[str]:
    """上传结构"""
    path = Path(file_path)
    content = path.read_text()
    
    response = requests.post(
        f"{BASE_URL}/structures/upload",
        data={
            "name": path.stem,
            "format": "cif",
            "content": content
        }
    )
    
    if response.ok:
        return response.json()['id']
    return None


def run_optimization_benchmark(
    structure_id: str,
    model_key: str,
    fmax: float = 0.05,
    max_steps: int = 100,
    timeout: int = 300
) -> Dict:
    """运行优化基准测试"""
    start_time = time.time()
    
    # 提交任务
    response = requests.post(
        f"{BASE_URL}/tasks/optimization",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "fmax": fmax,
            "max_steps": max_steps,
            "optimizer": "BFGS",
            "use_filter": True,
            "enable_d3": True
        }
    )
    
    if not response.ok:
        return {
            "model": model_key,
            "status": "submit_failed",
            "error": response.text
        }
    
    task_id = response.json()['task_id']
    
    # 等待完成
    while time.time() - start_time < timeout:
        resp = requests.get(f"{BASE_URL}/tasks/{task_id}")
        if not resp.ok:
            return {"model": model_key, "status": "error"}
        
        task = resp.json()
        
        if task['status'] == 'completed':
            result = requests.get(f"{BASE_URL}/results/{task_id}").json()
            wall_time = time.time() - start_time
            
            return {
                "model": model_key,
                "status": "completed",
                "final_energy": result.get('final_energy'),
                "energy_per_atom": result.get('energy_per_atom'),
                "num_steps": result.get('num_steps'),
                "converged": result.get('converged'),
                "rmsd": result.get('rmsd'),
                "wall_time": wall_time
            }
        
        elif task['status'] in ['failed', 'cancelled']:
            return {
                "model": model_key,
                "status": task['status'],
                "error": task.get('error_message')
            }
        
        time.sleep(2)
    
    return {"model": model_key, "status": "timeout"}


def run_full_benchmark(
    structure_paths: List[str],
    models: List[str] = None,
    output_file: str = None
):
    """
    运行完整基准测试
    
    Args:
        structure_paths: 结构文件路径列表
        models: 要测试的模型列表
        output_file: 结果输出文件
    """
    if models is None:
        models = QUICK_TEST_MODELS
    
    if output_file is None:
        output_file = f"benchmark_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    print(f"\n{'='*70}")
    print(f"MIRA 全模型基准测试")
    print(f"{'='*70}")
    print(f"测试结构: {len(structure_paths)} 个")
    print(f"测试模型: {len(models)} 个")
    print(f"输出文件: {output_file}")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    all_results = []
    
    for struct_path in structure_paths:
        struct_name = Path(struct_path).stem
        print(f"\n{'='*70}")
        print(f"结构: {struct_name}")
        print(f"{'='*70}")
        
        # 上传结构
        structure_id = upload_structure(struct_path)
        if not structure_id:
            print(f"  ✗ 上传失败")
            continue
        
        structure_results = {
            "structure": struct_name,
            "models": []
        }
        
        for model in models:
            print(f"\n  [{model}]", end=" ", flush=True)
            
            result = run_optimization_benchmark(
                structure_id=structure_id,
                model_key=model,
                fmax=0.05,
                max_steps=100
            )
            
            if result['status'] == 'completed':
                print(f"✓ E={result['energy_per_atom']:.4f} eV/atom, "
                      f"步数={result['num_steps']}, "
                      f"时间={result['wall_time']:.1f}s")
            else:
                print(f"✗ {result['status']}")
            
            structure_results['models'].append(result)
        
        all_results.append(structure_results)
    
    # 保存结果
    output_path = Path(__file__).parent / output_file
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n结果已保存到: {output_path}")
    
    # 生成汇总报告
    print_benchmark_summary(all_results)
    
    return all_results


def print_benchmark_summary(results: List[Dict]):
    """打印基准测试汇总"""
    print(f"\n{'='*70}")
    print("基准测试汇总")
    print(f"{'='*70}")
    
    # 收集所有模型的统计
    model_stats = {}
    
    for struct_result in results:
        for model_result in struct_result['models']:
            model = model_result['model']
            if model not in model_stats:
                model_stats[model] = {
                    'success': 0,
                    'failed': 0,
                    'energies': [],
                    'times': [],
                    'steps': []
                }
            
            if model_result['status'] == 'completed':
                model_stats[model]['success'] += 1
                model_stats[model]['energies'].append(model_result['energy_per_atom'])
                model_stats[model]['times'].append(model_result['wall_time'])
                model_stats[model]['steps'].append(model_result['num_steps'])
            else:
                model_stats[model]['failed'] += 1
    
    # 打印统计
    print(f"\n{'模型':<25} {'成功':<8} {'失败':<8} {'平均能量':<15} {'平均时间':<12} {'平均步数':<10}")
    print("-" * 80)
    
    for model, stats in model_stats.items():
        if stats['success'] > 0:
            avg_energy = sum(stats['energies']) / len(stats['energies'])
            avg_time = sum(stats['times']) / len(stats['times'])
            avg_steps = sum(stats['steps']) / len(stats['steps'])
            print(f"{model:<25} {stats['success']:<8} {stats['failed']:<8} "
                  f"{avg_energy:<15.4f} {avg_time:<12.1f}s {avg_steps:<10.1f}")
        else:
            print(f"{model:<25} {stats['success']:<8} {stats['failed']:<8} {'N/A':<15} {'N/A':<12} {'N/A':<10}")


def compare_energy_predictions(
    structure_path: str,
    models: List[str] = None,
    reference_energy: float = None
):
    """
    比较不同模型的能量预测
    
    Args:
        structure_path: 结构文件
        models: 模型列表
        reference_energy: 参考能量 (如 DFT 值)
    """
    if models is None:
        models = QUICK_TEST_MODELS
    
    print(f"\n{'='*60}")
    print(f"能量预测比较")
    print(f"{'='*60}")
    print(f"结构: {Path(structure_path).stem}")
    
    if reference_energy is not None:
        print(f"参考值: {reference_energy:.4f} eV/atom")
    
    structure_id = upload_structure(structure_path)
    if not structure_id:
        print("上传失败")
        return None
    
    results = []
    
    for model in models:
        print(f"\n  测试 {model}...", end=" ", flush=True)
        
        result = run_optimization_benchmark(
            structure_id=structure_id,
            model_key=model,
            fmax=0.01,
            max_steps=200
        )
        
        if result['status'] == 'completed':
            energy = result['energy_per_atom']
            error = None
            if reference_energy is not None:
                error = energy - reference_energy
            
            results.append({
                "model": model,
                "energy": energy,
                "error": error
            })
            
            if error is not None:
                print(f"E = {energy:.4f} eV/atom (误差: {error:+.4f})")
            else:
                print(f"E = {energy:.4f} eV/atom")
        else:
            results.append({
                "model": model,
                "energy": None,
                "error": None
            })
            print(f"失败: {result['status']}")
    
    # 汇总
    print(f"\n{'='*60}")
    print("能量预测汇总")
    print(f"{'='*60}")
    
    if reference_energy is not None:
        print(f"{'模型':<25} {'能量(eV/atom)':<18} {'误差':<15}")
        print("-" * 60)
        for r in results:
            if r['energy'] is not None:
                print(f"{r['model']:<25} {r['energy']:<18.4f} {r['error']:+.4f}")
    else:
        print(f"{'模型':<25} {'能量(eV/atom)':<18}")
        print("-" * 45)
        for r in results:
            if r['energy'] is not None:
                print(f"{r['model']:<25} {r['energy']:<18.4f}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    structures_dir = Path(__file__).parent / "structures"
    
    # 示例 1: 快速基准测试
    run_full_benchmark(
        structure_paths=[
            str(structures_dir / "HKUST-1.cif"),
            str(structures_dir / "ZIF-8.cif"),
            str(structures_dir / "UiO-66.cif"),
            str(structures_dir / "MOF-5.cif"),
        ],
        models=QUICK_TEST_MODELS,
        output_file="quick_benchmark.json"
    )
    
    # 示例 2: 完整基准测试（所有模型）
    # 警告：这可能需要很长时间！
    # run_full_benchmark(
    #     structure_paths=[...],
    #     models=ALL_MODELS,
    #     output_file="full_benchmark.json"
    # )
    
    # 示例 3: 能量预测比较
    # compare_energy_predictions(
    #     str(structures_dir / "HKUST-1.cif"),
    #     models=QUICK_TEST_MODELS,
    #     reference_energy=-5.234  # 假设的 DFT 参考值
    # )
