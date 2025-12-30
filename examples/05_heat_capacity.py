"""
MIRA Examples - 热容计算
使用 Phonopy 进行声子计算，获得热容和其他热力学性质

运行前确保:
1. 已安装 ML 力场: python scripts/install_models.py --check
2. 服务已启动: uvicorn app.main:app --host 0.0.0.0 --port 8000
"""
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
BASE_URL = "http://localhost:8000/api/v1"


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


def submit_heat_capacity(
    structure_id: str,
    model_key: str,
    supercell: list = None,
    mesh: list = None,
    temperature_range: tuple = (0, 1000),
    temperature_points: int = 101,
    pre_optimize: bool = True,
    enable_d3: bool = True
) -> Optional[str]:
    """
    提交热容计算任务
    
    Args:
        structure_id: 结构 ID
        model_key: 模型键
        supercell: 超胞大小 [nx, ny, nz]
        mesh: 声子 q 点网格 [qx, qy, qz]
        temperature_range: 温度范围 (Tmin, Tmax) in K
        temperature_points: 温度点数
        pre_optimize: 是否预优化
        enable_d3: 启用 D3 校正
        
    Returns:
        task_id
    """
    if supercell is None:
        supercell = [1, 1, 1]  # MOF 晶胞较大，通常不需要大超胞
    if mesh is None:
        mesh = [2, 2, 2]  # 较稀疏的 mesh，加快计算
    
    response = requests.post(
        f"{BASE_URL}/tasks/heat-capacity",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "supercell": supercell,
            "mesh": mesh,
            "temperature_range": temperature_range,
            "temperature_points": temperature_points,
            "displacement": 0.01,
            "pre_optimize": pre_optimize,
            "enable_d3": enable_d3
        }
    )
    
    if response.ok:
        return response.json()['task_id']
    else:
        print(f"提交失败: {response.text}")
        return None


def wait_for_task(task_id: str, timeout: int = 1800) -> dict:
    """等待任务完成（声子计算可能较慢）"""
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
        print(f"\r  声子计算: {progress:.1f}% (已用时 {elapsed:.0f}s)", end="", flush=True)
        
        time.sleep(5)
    
    return {"status": "timeout"}


def calculate_heat_capacity(
    structure_path: str,
    model_key: str = "mace-mp",
    supercell: list = None,
    mesh: list = None
):
    """
    计算 MOF 的热容
    
    Args:
        structure_path: 结构文件路径
        model_key: 模型
        supercell: 超胞大小
        mesh: q 点网格
    """
    if supercell is None:
        supercell = [1, 1, 1]
    if mesh is None:
        mesh = [2, 2, 2]
    
    print(f"\n{'='*60}")
    print(f"热容计算 (声子方法)")
    print(f"{'='*60}")
    print(f"结构: {Path(structure_path).stem}")
    print(f"模型: {model_key}")
    print(f"超胞: {supercell}")
    print(f"q 点网格: {mesh}")
    
    # 上传结构
    print("\n上传结构...")
    structure_id = upload_structure(structure_path)
    if not structure_id:
        print("上传失败")
        return None
    
    # 提交任务
    print("提交声子计算任务...")
    task_id = submit_heat_capacity(
        structure_id=structure_id,
        model_key=model_key,
        supercell=supercell,
        mesh=mesh
    )
    
    if not task_id:
        return None
    
    print(f"任务 ID: {task_id}")
    print("等待计算完成（声子计算可能需要较长时间）...")
    
    # 等待结果
    result = wait_for_task(task_id, timeout=3600)
    print()
    
    if result.get('status') in ['failed', 'timeout', 'error']:
        print(f"计算失败: {result.get('error', result.get('status'))}")
        return None
    
    # 显示结果
    print(f"\n{'='*60}")
    print("热容计算结果")
    print(f"{'='*60}")
    
    cv_300k = result.get('cv_300k', 0)
    entropy_300k = result.get('entropy_300k', 0)
    free_energy_300k = result.get('free_energy_300k', 0)
    internal_energy_300k = result.get('internal_energy_300k', 0)
    
    print(f"\n300K 热力学性质:")
    print(f"  热容 Cv: {cv_300k:.4f} J/(mol·K)")
    print(f"  熵 S: {entropy_300k:.4f} J/(mol·K)")
    print(f"  亥姆霍兹自由能 F: {free_energy_300k:.4f} kJ/mol")
    print(f"  内能 U: {internal_energy_300k:.4f} kJ/mol")
    
    # 虚频检查
    has_imaginary = result.get('has_imaginary_modes', False)
    num_imaginary = result.get('num_imaginary_modes', 0)
    
    print(f"\n声子稳定性:")
    if has_imaginary:
        print(f"  ⚠ 存在 {num_imaginary} 个虚频模式")
        print(f"    这可能表示结构未完全优化或动力学不稳定")
    else:
        print(f"  ✓ 无虚频，声子稳定")
    
    # 显示部分 T-Cv 数据
    temp_cv_data = result.get('temperature_cv_data', [])
    if temp_cv_data:
        print(f"\n温度-热容曲线 (部分数据):")
        print(f"{'温度(K)':<12} {'Cv (J/mol·K)':<15}")
        print("-" * 30)
        
        # 显示特定温度点
        target_temps = [100, 200, 300, 400, 500, 600, 800]
        for point in temp_cv_data:
            T = point['temperature']
            if any(abs(T - t) < 10 for t in target_temps):
                print(f"{T:<12.0f} {point['cv']:<15.4f}")
    
    # 显示部分声子频率
    frequencies = result.get('phonon_frequencies', [])
    if frequencies:
        print(f"\n声子频率 (前 20 个, THz):")
        for i, freq in enumerate(frequencies[:20]):
            marker = " ⚠" if freq < 0 else ""
            print(f"  {i+1:2d}: {freq:8.4f}{marker}")
    
    return result


def compare_cv_across_mofs(
    structure_paths: list,
    model_key: str = "mace-mp"
):
    """
    比较不同 MOF 的热容
    
    Args:
        structure_paths: 结构文件路径列表
        model_key: 模型
    """
    print(f"\n{'='*60}")
    print(f"MOF 热容比较")
    print(f"模型: {model_key}")
    print(f"{'='*60}")
    
    results = []
    
    for path in structure_paths:
        name = Path(path).stem
        print(f"\n--- {name} ---")
        
        structure_id = upload_structure(path)
        if not structure_id:
            results.append({"name": name, "status": "upload_failed"})
            continue
        
        task_id = submit_heat_capacity(
            structure_id=structure_id,
            model_key=model_key,
            supercell=[1, 1, 1],
            mesh=[2, 2, 2]
        )
        
        if not task_id:
            results.append({"name": name, "status": "submit_failed"})
            continue
        
        result = wait_for_task(task_id)
        print()
        
        if result.get('status') in ['failed', 'timeout', 'error']:
            results.append({"name": name, "status": result.get('status')})
        else:
            results.append({
                "name": name,
                "status": "completed",
                "cv_300k": result.get('cv_300k', 0),
                "entropy_300k": result.get('entropy_300k', 0),
                "has_imaginary": result.get('has_imaginary_modes', False)
            })
            print(f"  Cv(300K) = {result.get('cv_300k', 0):.4f} J/(mol·K)")
    
    # 汇总
    print(f"\n{'='*60}")
    print("热容比较汇总")
    print(f"{'='*60}")
    print(f"{'MOF':<15} {'Cv(300K)':<15} {'S(300K)':<15} {'虚频':<8}")
    print("-" * 55)
    
    for r in results:
        if r['status'] == 'completed':
            imaginary = "是" if r['has_imaginary'] else "否"
            print(f"{r['name']:<15} {r['cv_300k']:<15.4f} {r['entropy_300k']:<15.4f} {imaginary:<8}")
        else:
            print(f"{r['name']:<15} {'失败'}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    structures_dir = Path(__file__).parent / "structures"
    
    # 示例 1: 计算 ZIF-8 的热容
    # ZIF-8 晶胞较小，适合声子计算
    calculate_heat_capacity(
        str(structures_dir / "ZIF-8.cif"),
        model_key="mace-mp",
        supercell=[1, 1, 1],
        mesh=[2, 2, 2]
    )
    
    # 示例 2: 比较不同 MOF 的热容
    # compare_cv_across_mofs(
    #     [
    #         str(structures_dir / "ZIF-8.cif"),
    #         str(structures_dir / "UiO-66.cif"),
    #     ],
    #     model_key="mace-mp"
    # )
