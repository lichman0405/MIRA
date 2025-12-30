"""
MIRA Examples - 体积模量计算
通过 E-V 曲线拟合 Birch-Murnaghan 状态方程计算体积模量

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


def submit_bulk_modulus(
    structure_id: str,
    model_key: str,
    volume_range: float = 0.1,
    num_points: int = 11,
    pre_optimize: bool = True,
    enable_d3: bool = True
) -> Optional[str]:
    """
    提交体积模量计算任务
    
    Args:
        structure_id: 结构 ID
        model_key: 模型键
        volume_range: 体积变化范围 (±比例)
        num_points: E-V 曲线采样点数
        pre_optimize: 是否预优化
        enable_d3: 启用 D3 校正
        
    Returns:
        task_id
    """
    response = requests.post(
        f"{BASE_URL}/tasks/bulk-modulus",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "volume_range": volume_range,
            "num_points": num_points,
            "pre_optimize": pre_optimize,
            "enable_d3": enable_d3
        }
    )
    
    if response.ok:
        return response.json()['task_id']
    else:
        print(f"提交失败: {response.text}")
        return None


def wait_for_task(task_id: str, timeout: int = 600) -> dict:
    """等待任务完成"""
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
        print(f"\r  E-V 曲线计算: {progress:.1f}%", end="", flush=True)
        
        time.sleep(2)
    
    return {"status": "timeout"}


def calculate_bulk_modulus(
    structure_path: str,
    model_key: str = "mace-mp",
    volume_range: float = 0.08,
    num_points: int = 11
):
    """
    计算 MOF 的体积模量
    
    Args:
        structure_path: 结构文件路径
        model_key: 模型
        volume_range: 体积变化范围
        num_points: 采样点数
    """
    print(f"\n{'='*60}")
    print(f"体积模量计算")
    print(f"{'='*60}")
    print(f"结构: {Path(structure_path).stem}")
    print(f"模型: {model_key}")
    print(f"体积范围: ±{volume_range*100}%")
    print(f"采样点数: {num_points}")
    
    # 上传结构
    print("\n上传结构...")
    structure_id = upload_structure(structure_path)
    if not structure_id:
        print("上传失败")
        return None
    
    # 提交任务
    print("提交计算任务...")
    task_id = submit_bulk_modulus(
        structure_id=structure_id,
        model_key=model_key,
        volume_range=volume_range,
        num_points=num_points
    )
    
    if not task_id:
        return None
    
    print(f"任务 ID: {task_id}")
    
    # 等待结果
    result = wait_for_task(task_id)
    print()
    
    if result.get('status') in ['failed', 'timeout', 'error']:
        print(f"计算失败: {result.get('error', result.get('status'))}")
        return None
    
    # 显示结果
    print(f"\n{'='*60}")
    print("体积模量计算结果")
    print(f"{'='*60}")
    
    B0 = result.get('bulk_modulus', 0)
    B0_prime = result.get('bulk_modulus_derivative', 0)
    V0 = result.get('equilibrium_volume', 0)
    E0 = result.get('equilibrium_energy', 0)
    r2 = result.get('r_squared', 0)
    
    print(f"\n体积模量 (B₀): {B0:.2f} GPa")
    print(f"B₀ 导数 (B₀'): {B0_prime:.2f}")
    print(f"平衡体积 (V₀): {V0:.2f} Å³")
    print(f"平衡能量 (E₀): {E0:.4f} eV")
    print(f"拟合质量 (R²): {r2:.6f}")
    
    # E-V 曲线数据
    ev_data = result.get('ev_curve_data', [])
    if ev_data:
        print(f"\nE-V 曲线数据 ({len(ev_data)} 点):")
        print(f"{'体积(Å³)':<15} {'能量(eV)':<15}")
        print("-" * 30)
        for point in ev_data:
            print(f"{point['volume']:<15.2f} {point['energy']:<15.4f}")
    
    # 力学稳定性评估
    print(f"\n力学稳定性评估:")
    if B0 > 10:
        print(f"  ✓ 体积模量较高 ({B0:.1f} GPa)，力学稳定性好")
    elif B0 > 2:
        print(f"  △ 体积模量中等 ({B0:.1f} GPa)，力学稳定性一般")
    else:
        print(f"  ✗ 体积模量较低 ({B0:.1f} GPa)，力学稳定性可能较差")
    
    return result


def compare_models_bulk_modulus(
    structure_path: str,
    models: list = None
):
    """
    使用不同模型计算体积模量并比较
    
    Args:
        structure_path: 结构文件路径
        models: 模型列表
    """
    if models is None:
        models = ["mace-mp", "orb-v2", "grace-2l", "sevennet-0"]
    
    print(f"\n{'='*60}")
    print(f"多模型体积模量比较")
    print(f"结构: {Path(structure_path).stem}")
    print(f"{'='*60}")
    
    results = []
    
    for model in models:
        print(f"\n--- 模型: {model} ---")
        
        # 上传结构（每个模型单独上传避免冲突）
        structure_id = upload_structure(structure_path)
        if not structure_id:
            continue
        
        # 提交任务
        task_id = submit_bulk_modulus(
            structure_id=structure_id,
            model_key=model
        )
        
        if not task_id:
            results.append({"model": model, "status": "failed"})
            continue
        
        # 等待结果
        result = wait_for_task(task_id)
        print()
        
        if result.get('status') in ['failed', 'timeout', 'error']:
            results.append({"model": model, "status": result.get('status')})
        else:
            results.append({
                "model": model,
                "status": "completed",
                "bulk_modulus": result.get('bulk_modulus', 0),
                "r_squared": result.get('r_squared', 0)
            })
            print(f"  B₀ = {result.get('bulk_modulus', 0):.2f} GPa (R² = {result.get('r_squared', 0):.4f})")
    
    # 汇总比较
    print(f"\n{'='*60}")
    print("体积模量比较汇总")
    print(f"{'='*60}")
    print(f"{'模型':<20} {'B₀ (GPa)':<15} {'R²':<10}")
    print("-" * 45)
    
    for r in results:
        if r['status'] == 'completed':
            print(f"{r['model']:<20} {r['bulk_modulus']:<15.2f} {r['r_squared']:<10.6f}")
        else:
            print(f"{r['model']:<20} {'失败'}")
    
    return results


# ========== 主程序 ==========
if __name__ == "__main__":
    structures_dir = Path(__file__).parent / "structures"
    
    # 示例 1: 单模型计算 UiO-66 的体积模量
    # UiO-66 是一个力学稳定性很好的 MOF
    calculate_bulk_modulus(
        str(structures_dir / "UiO-66.cif"),
        model_key="mace-mp",
        volume_range=0.08,
        num_points=11
    )
    
    # 示例 2: 多模型比较
    # compare_models_bulk_modulus(
    #     str(structures_dir / "UiO-66.cif"),
    #     models=["mace-mp", "orb-v2", "grace-2l"]
    # )
