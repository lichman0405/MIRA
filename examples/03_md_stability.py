"""
MIRA Examples - MD 稳定性测试
对 MOF 结构进行分子动力学模拟，评估热稳定性

运行前确保:
1. 已安装 ML 力场: python scripts/install_models.py --check
2. 服务已启动: 
   - 本地: uvicorn app.main:app --host 0.0.0.0 --port 8000
   - Docker: ./scripts/deploy.sh test

支持环境变量:
    MIRA_GATEWAY_URL=http://192.168.100.207:8000  # 测试服务器
"""
import os
import requests
import time
from pathlib import Path
from typing import Optional

# 依赖检查
try:
    from setup_check import ensure_dependencies, get_first_available_model
    ensure_dependencies(verbose=False)
except ImportError:
    print("提示: 运行 'python scripts/install_models.py --check' 检查依赖")

# ========== 配置 ==========
GATEWAY_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000")
BASE_URL = f"{GATEWAY_URL}/api/v1"


def upload_structure(file_path: str) -> Optional[str]:
    """上传结构文件"""
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


def submit_stability_test(
    structure_id: str,
    model_key: str,
    nvt_steps: int = 2000,
    nvt_temperature: float = 300.0,
    npt_steps: int = 5000,
    npt_temperature: float = 300.0,
    npt_pressure: float = 1.01325,
    timestep: float = 1.0,
    pre_optimize: bool = True,
    enable_d3: bool = True
) -> Optional[str]:
    """
    提交 MD 稳定性测试任务
    
    Args:
        structure_id: 结构 ID
        model_key: 模型键
        nvt_steps: NVT 平衡步数
        nvt_temperature: NVT 温度 (K)
        npt_steps: NPT 生产运行步数
        npt_temperature: NPT 温度 (K)
        npt_pressure: 压力 (bar)
        timestep: 时间步长 (fs)
        pre_optimize: 是否预优化
        enable_d3: 启用 D3 校正
        
    Returns:
        task_id
    """
    response = requests.post(
        f"{BASE_URL}/tasks/stability",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "nvt_steps": nvt_steps,
            "nvt_temperature": nvt_temperature,
            "npt_steps": npt_steps,
            "npt_temperature": npt_temperature,
            "npt_pressure": npt_pressure,
            "timestep": timestep,
            "pre_optimize": pre_optimize,
            "enable_d3": enable_d3,
            "trajectory_interval": 50
        }
    )
    
    if response.ok:
        return response.json()['task_id']
    else:
        print(f"提交失败: {response.text}")
        return None


def wait_for_task(task_id: str, timeout: int = 1800) -> dict:
    """等待任务完成（MD 任务可能需要较长时间）"""
    start_time = time.time()
    
    while time.time() - start_time < timeout:
        response = requests.get(f"{BASE_URL}/tasks/{task_id}")
        if not response.ok:
            return {"status": "error"}
        
        task = response.json()
        status = task['status']
        
        if status == "completed":
            result_response = requests.get(f"{BASE_URL}/results/{task_id}")
            return result_response.json()
        
        elif status in ["failed", "cancelled"]:
            return {"status": status, "error": task.get('error_message')}
        
        progress = task.get('progress', 0)
        elapsed = time.time() - start_time
        print(f"\r  进度: {progress:.1f}% (已用时 {elapsed:.0f}s)", end="", flush=True)
        
        time.sleep(5)
    
    return {"status": "timeout"}


def run_stability_test(
    structure_path: str,
    model_key: str = "mace-mp",
    temperature: float = 300.0,
    total_time_ps: float = 5.0
):
    """
    运行 MD 稳定性测试
    
    Args:
        structure_path: 结构文件路径
        model_key: 使用的模型
        temperature: 目标温度 (K)
        total_time_ps: 总模拟时间 (ps)
    """
    timestep = 1.0  # fs
    total_steps = int(total_time_ps * 1000 / timestep)
    nvt_steps = total_steps // 4
    npt_steps = total_steps - nvt_steps
    
    print(f"\n{'='*60}")
    print(f"MD 稳定性测试")
    print(f"{'='*60}")
    print(f"结构: {Path(structure_path).stem}")
    print(f"模型: {model_key}")
    print(f"温度: {temperature} K")
    print(f"总时间: {total_time_ps} ps")
    print(f"  NVT 步数: {nvt_steps}")
    print(f"  NPT 步数: {npt_steps}")
    
    # 上传结构
    print("\n上传结构...")
    structure_id = upload_structure(structure_path)
    if not structure_id:
        print("上传失败")
        return None
    
    # 提交任务
    print("提交 MD 任务...")
    task_id = submit_stability_test(
        structure_id=structure_id,
        model_key=model_key,
        nvt_steps=nvt_steps,
        nvt_temperature=temperature,
        npt_steps=npt_steps,
        npt_temperature=temperature
    )
    
    if not task_id:
        return None
    
    print(f"任务 ID: {task_id}")
    print("等待任务完成...")
    
    # 等待结果
    result = wait_for_task(task_id, timeout=3600)
    print()
    
    if result.get('status') in ['failed', 'timeout', 'error']:
        print(f"任务失败: {result.get('error', result.get('status'))}")
        return None
    
    # 显示结果
    print(f"\n{'='*60}")
    print("稳定性测试结果")
    print(f"{'='*60}")
    print(f"模拟时间: {result.get('simulation_time_ps', 0):.2f} ps")
    print(f"总步数: {result.get('total_steps', 0)}")
    print(f"\n温度统计:")
    print(f"  平均温度: {result.get('average_temperature', 0):.1f} K")
    print(f"  温度标准差: {result.get('temperature_std', 0):.1f} K")
    print(f"\n结构变化:")
    print(f"  最终 RMSD: {result.get('final_rmsd', 0):.4f} Å")
    print(f"  最大 RMSD: {result.get('max_rmsd', 0):.4f} Å")
    print(f"  体积漂移: {result.get('volume_drift_percent', 0):.2f}%")
    print(f"\n能量统计:")
    print(f"  平均能量: {result.get('average_energy', 0):.4f} eV")
    print(f"  能量标准差: {result.get('energy_std', 0):.4f} eV")
    
    # 配位数分析
    coord = result.get('coordination_analysis', {})
    if coord:
        print(f"\n配位数分析:")
        initial = coord.get('initial', {})
        final = coord.get('final', {})
        for metal in initial:
            init_avg = sum(initial[metal]) / len(initial[metal]) if initial[metal] else 0
            final_avg = sum(final.get(metal, [])) / len(final.get(metal, [])) if final.get(metal) else 0
            print(f"  {metal}: {init_avg:.1f} -> {final_avg:.1f}")
    
    # 稳定性评估
    print(f"\n稳定性评估:")
    final_rmsd = result.get('final_rmsd', 0)
    volume_drift = abs(result.get('volume_drift_percent', 0))
    
    if final_rmsd < 0.5 and volume_drift < 5:
        print("  ✓ 结构高度稳定")
    elif final_rmsd < 1.0 and volume_drift < 10:
        print("  △ 结构基本稳定")
    else:
        print("  ✗ 结构可能不稳定")
    
    return result


def run_temperature_series(
    structure_path: str,
    model_key: str = "mace-mp",
    temperatures: list = None
):
    """
    运行温度系列测试
    
    Args:
        structure_path: 结构文件路径
        model_key: 模型
        temperatures: 温度列表 (K)
    """
    if temperatures is None:
        temperatures = [300, 400, 500, 600]
    
    print(f"\n{'='*60}")
    print(f"温度系列稳定性测试")
    print(f"结构: {Path(structure_path).stem}")
    print(f"模型: {model_key}")
    print(f"温度: {temperatures} K")
    print(f"{'='*60}")
    
    results = []
    
    for T in temperatures:
        print(f"\n--- 温度: {T} K ---")
        result = run_stability_test(
            structure_path,
            model_key=model_key,
            temperature=T,
            total_time_ps=2.0  # 较短的模拟时间
        )
        
        if result:
            results.append({
                "temperature": T,
                "final_rmsd": result.get('final_rmsd', 0),
                "volume_drift": result.get('volume_drift_percent', 0)
            })
    
    # 汇总
    print(f"\n{'='*60}")
    print("温度系列结果汇总")
    print(f"{'='*60}")
    print(f"{'温度(K)':<12} {'RMSD(Å)':<12} {'体积漂移(%)':<12}")
    print("-" * 40)
    
    for r in results:
        print(f"{r['temperature']:<12} {r['final_rmsd']:<12.4f} {r['volume_drift']:<12.2f}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    structures_dir = Path(__file__).parent / "structures"
    
    # 示例 1: 单个结构的稳定性测试
    run_stability_test(
        str(structures_dir / "HKUST-1.cif"),
        model_key="mace-mp",
        temperature=300.0,
        total_time_ps=5.0
    )
    
    # 示例 2: 温度系列测试（评估热稳定性范围）
    # run_temperature_series(
    #     str(structures_dir / "ZIF-8.cif"),
    #     model_key="mace-mp",
    #     temperatures=[300, 400, 500]
    # )
