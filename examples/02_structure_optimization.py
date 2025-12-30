"""
MIRA Examples - 结构优化示例
使用不同 ML 力场模型对 MOF 结构进行优化

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
from pathlib import Path
from typing import Optional

# 客户端初始化
from client_utils import init_client, get_service_models

GATEWAY_URL, BASE_URL = init_client(verbose=True)

# 所有可用模型（按家族分组）
ALL_MODELS = {
    "MACE": ["mace-mp", "mace-off23", "mace-omat", "mace-mpa", "mace-ani"],
    "ORB": ["orb-v2", "orb-d3-v2", "orb-v3", "orb-v3-mpa", "orb-omat-v3-lora"],
    "OMAT24": ["omat24-base", "omat24-large", "eqv2-omat", "eqv2-mptrj"],
    "GRACE": ["grace-2l", "grace-2l-omat", "grace-2m"],
    "MatterSim": ["mattersim-5m"],
    "SevenNet": ["sevennet-0", "sevennet-mf-ompa", "sevennet-l3i5"],
    "PosEGNN": ["posegnn"],
    "MatGL": ["m3gnet", "chgnet"]
}


def upload_structure(file_path: str) -> Optional[str]:
    """上传结构文件，返回 structure_id"""
    path = Path(file_path)
    content = path.read_text()
    
    ext = path.suffix.lower()
    format_map = {'.cif': 'cif', '.xyz': 'xyz'}
    
    response = requests.post(
        f"{BASE_URL}/structures/upload",
        data={
            "name": path.stem,
            "format": format_map.get(ext, 'cif'),
            "content": content
        }
    )
    
    if response.ok:
        return response.json()['id']
    return None


def submit_optimization(
    structure_id: str,
    model_key: str,
    fmax: float = 0.05,
    max_steps: int = 200,
    optimizer: str = "BFGS",
    enable_d3: bool = True
) -> Optional[str]:
    """
    提交结构优化任务
    
    Args:
        structure_id: 结构 ID
        model_key: 模型键
        fmax: 力收敛阈值 (eV/Å)
        max_steps: 最大优化步数
        optimizer: 优化器 (BFGS/FIRE/LBFGS)
        enable_d3: 是否启用 D3 色散校正
        
    Returns:
        task_id
    """
    response = requests.post(
        f"{BASE_URL}/tasks/optimization",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "fmax": fmax,
            "max_steps": max_steps,
            "optimizer": optimizer,
            "use_filter": True,  # 使用 FrechetCellFilter
            "enable_d3": enable_d3
        }
    )
    
    if response.ok:
        return response.json()['task_id']
    else:
        print(f"提交失败 ({model_key}): {response.text}")
        return None


def wait_for_task(task_id: str, timeout: int = 600) -> dict:
    """等待任务完成"""
    start_time = time.time()
    
    while time.time() - start_time < timeout:
        response = requests.get(f"{BASE_URL}/tasks/{task_id}")
        if not response.ok:
            return {"status": "error", "error": response.text}
        
        task = response.json()
        status = task['status']
        
        if status == "completed":
            # 获取结果
            result_response = requests.get(f"{BASE_URL}/results/{task_id}")
            return result_response.json()
        
        elif status == "failed":
            return {"status": "failed", "error": task.get('error_message')}
        
        elif status == "cancelled":
            return {"status": "cancelled"}
        
        # 显示进度
        progress = task.get('progress', 0)
        print(f"\r  进度: {progress:.1f}%", end="", flush=True)
        
        time.sleep(2)
    
    return {"status": "timeout"}


def run_optimization_benchmark(
    structure_path: str,
    models: list = None,
    fmax: float = 0.05
):
    """
    对单个结构运行多模型优化基准测试
    
    Args:
        structure_path: 结构文件路径
        models: 要测试的模型列表（默认测试部分代表性模型）
        fmax: 力收敛阈值
    """
    # 默认选择每个家族的代表性模型
    if models is None:
        models = [
            "mace-mp",
            "orb-v2", 
            "omat24-base",
            "grace-2l",
            "mattersim-5m",
            "sevennet-0"
        ]
    
    print(f"\n{'='*60}")
    print(f"结构优化基准测试: {Path(structure_path).stem}")
    print(f"{'='*60}")
    
    # 上传结构
    print("\n上传结构...")
    structure_id = upload_structure(structure_path)
    if not structure_id:
        print("结构上传失败")
        return
    print(f"结构 ID: {structure_id}")
    
    # 存储结果
    results = []
    
    # 对每个模型运行优化
    for model in models:
        print(f"\n--- 模型: {model} ---")
        
        # 提交任务
        task_id = submit_optimization(
            structure_id=structure_id,
            model_key=model,
            fmax=fmax
        )
        
        if not task_id:
            results.append({
                "model": model,
                "status": "submit_failed"
            })
            continue
        
        print(f"任务 ID: {task_id}")
        
        # 等待完成
        result = wait_for_task(task_id)
        print()  # 换行
        
        if result.get('status') in ['failed', 'timeout', 'cancelled', 'error']:
            results.append({
                "model": model,
                "status": result.get('status'),
                "error": result.get('error')
            })
            print(f"  状态: {result.get('status')}")
        else:
            results.append({
                "model": model,
                "status": "completed",
                "final_energy": result.get('final_energy'),
                "energy_per_atom": result.get('energy_per_atom'),
                "num_steps": result.get('num_steps'),
                "converged": result.get('converged'),
                "rmsd": result.get('rmsd'),
                "volume_change": result.get('volume_change_percent')
            })
            print(f"  最终能量: {result.get('final_energy', 0):.4f} eV")
            print(f"  每原子能量: {result.get('energy_per_atom', 0):.4f} eV/atom")
            print(f"  优化步数: {result.get('num_steps')}")
            print(f"  收敛: {result.get('converged')}")
            print(f"  RMSD: {result.get('rmsd', 0):.4f} Å")
            print(f"  体积变化: {result.get('volume_change_percent', 0):.2f}%")
    
    # 汇总结果
    print(f"\n{'='*60}")
    print("优化结果汇总")
    print(f"{'='*60}")
    print(f"{'模型':<20} {'能量(eV/atom)':<15} {'步数':<8} {'收敛':<6} {'RMSD(Å)':<10}")
    print("-" * 60)
    
    for r in results:
        if r['status'] == 'completed':
            print(f"{r['model']:<20} {r['energy_per_atom']:<15.4f} {r['num_steps']:<8} "
                  f"{'是' if r['converged'] else '否':<6} {r['rmsd']:<10.4f}")
        else:
            print(f"{r['model']:<20} {'失败 - ' + r['status']}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    # 示例结构目录
    structures_dir = Path(__file__).parent / "structures"
    
    # 选择要测试的 MOF 结构
    test_structures = [
        structures_dir / "HKUST-1.cif",
        structures_dir / "ZIF-8.cif",
    ]
    
    # 选择要测试的模型子集
    test_models = [
        "mace-mp",
        "mace-omat",
        "orb-v2",
        "grace-2l",
    ]
    
    for structure_path in test_structures:
        if structure_path.exists():
            run_optimization_benchmark(
                str(structure_path),
                models=test_models,
                fmax=0.05
            )
        else:
            print(f"结构文件不存在: {structure_path}")
