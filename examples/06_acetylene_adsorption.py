"""
MIRA Examples - 乙炔吸附相互作用能分析
计算 MOF 与乙炔 (C2H2) 分子之间的相互作用能

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
from client_utils import init_client

GATEWAY_URL, BASE_URL = init_client(verbose=True)


def upload_structure(file_path: str, name: str = None) -> Optional[str]:
    """上传结构文件"""
    path = Path(file_path)
    content = path.read_text()
    
    ext = path.suffix.lower()
    format_map = {'.cif': 'cif', '.xyz': 'xyz'}
    
    response = requests.post(
        f"{BASE_URL}/structures/upload",
        data={
            "name": name or path.stem,
            "format": format_map.get(ext, 'cif'),
            "content": content
        }
    )
    
    if response.ok:
        return response.json()['id']
    return None


def calculate_single_point_energy(
    structure_id: str,
    model_key: str,
    enable_d3: bool = True
) -> Optional[float]:
    """
    计算单点能量
    
    使用优化任务但 max_steps=0 来获取单点能量
    """
    response = requests.post(
        f"{BASE_URL}/tasks/optimization",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "fmax": 0.01,
            "max_steps": 1,  # 只计算初始能量
            "optimizer": "BFGS",
            "use_filter": False,
            "enable_d3": enable_d3
        }
    )
    
    if not response.ok:
        return None
    
    task_id = response.json()['task_id']
    
    # 等待
    for _ in range(60):
        resp = requests.get(f"{BASE_URL}/tasks/{task_id}")
        if not resp.ok:
            return None
        
        task = resp.json()
        if task['status'] == 'completed':
            result = requests.get(f"{BASE_URL}/results/{task_id}").json()
            return result.get('initial_energy')
        elif task['status'] in ['failed', 'cancelled']:
            return None
        
        time.sleep(1)
    
    return None


def calculate_interaction_energy_manual(
    host_structure_path: str,
    guest_structure_path: str,
    complex_structure_path: str,
    model_key: str = "mace-mp",
    enable_d3: bool = True
):
    """
    手动计算相互作用能
    
    E_interaction = E_complex - E_host - E_guest
    
    Args:
        host_structure_path: 宿主 (MOF) 结构文件
        guest_structure_path: 客体 (乙炔) 结构文件
        complex_structure_path: 复合物 (MOF+乙炔) 结构文件
        model_key: 模型
        enable_d3: 启用 D3 色散校正
    """
    print(f"\n{'='*60}")
    print(f"乙炔吸附相互作用能计算")
    print(f"{'='*60}")
    print(f"模型: {model_key}")
    print(f"D3 校正: {'是' if enable_d3 else '否'}")
    print(f"\n组分:")
    print(f"  宿主 (MOF): {Path(host_structure_path).stem}")
    print(f"  客体 (C2H2): {Path(guest_structure_path).stem}")
    print(f"  复合物: {Path(complex_structure_path).stem}")
    
    # 上传结构
    print("\n上传结构...")
    host_id = upload_structure(host_structure_path, "host_MOF")
    guest_id = upload_structure(guest_structure_path, "guest_C2H2")
    complex_id = upload_structure(complex_structure_path, "complex_MOF_C2H2")
    
    if not all([host_id, guest_id, complex_id]):
        print("结构上传失败")
        return None
    
    print("结构上传成功")
    
    # 计算各组分能量
    print("\n计算能量...")
    
    print("  计算复合物能量...", end=" ", flush=True)
    E_complex = calculate_single_point_energy(complex_id, model_key, enable_d3)
    if E_complex is None:
        print("失败")
        return None
    print(f"{E_complex:.4f} eV")
    
    print("  计算宿主能量...", end=" ", flush=True)
    E_host = calculate_single_point_energy(host_id, model_key, enable_d3)
    if E_host is None:
        print("失败")
        return None
    print(f"{E_host:.4f} eV")
    
    print("  计算客体能量...", end=" ", flush=True)
    E_guest = calculate_single_point_energy(guest_id, model_key, enable_d3)
    if E_guest is None:
        print("失败")
        return None
    print(f"{E_guest:.4f} eV")
    
    # 计算相互作用能
    E_interaction = E_complex - E_host - E_guest
    
    # 转换单位
    E_interaction_kJ_mol = E_interaction * 96.485  # eV -> kJ/mol
    E_interaction_kcal_mol = E_interaction * 23.061  # eV -> kcal/mol
    
    # 显示结果
    print(f"\n{'='*60}")
    print("相互作用能计算结果")
    print(f"{'='*60}")
    print(f"\n能量分解:")
    print(f"  E(复合物) = {E_complex:.4f} eV")
    print(f"  E(宿主)   = {E_host:.4f} eV")
    print(f"  E(客体)   = {E_guest:.4f} eV")
    print(f"\n相互作用能 = E(复合物) - E(宿主) - E(客体)")
    print(f"\n  ΔE = {E_interaction:.4f} eV")
    print(f"     = {E_interaction_kJ_mol:.2f} kJ/mol")
    print(f"     = {E_interaction_kcal_mol:.2f} kcal/mol")
    
    # 结合强度评估
    print(f"\n吸附强度评估:")
    if E_interaction_kJ_mol < -50:
        print(f"  ✓ 强吸附 (ΔE < -50 kJ/mol)")
    elif E_interaction_kJ_mol < -20:
        print(f"  △ 中等吸附 (-50 < ΔE < -20 kJ/mol)")
    elif E_interaction_kJ_mol < 0:
        print(f"  ○ 弱吸附 (-20 < ΔE < 0 kJ/mol)")
    else:
        print(f"  ✗ 排斥作用 (ΔE > 0)")
    
    return {
        "E_complex": E_complex,
        "E_host": E_host,
        "E_guest": E_guest,
        "E_interaction_eV": E_interaction,
        "E_interaction_kJ_mol": E_interaction_kJ_mol,
        "E_interaction_kcal_mol": E_interaction_kcal_mol
    }


def compare_models_interaction_energy(
    host_path: str,
    guest_path: str,
    complex_path: str,
    models: list = None
):
    """
    使用不同模型比较相互作用能
    """
    if models is None:
        models = ["mace-mp", "mace-omat", "orb-v2", "grace-2l"]
    
    print(f"\n{'='*60}")
    print(f"多模型相互作用能比较")
    print(f"{'='*60}")
    
    results = []
    
    for model in models:
        print(f"\n--- 模型: {model} ---")
        
        result = calculate_interaction_energy_manual(
            host_path, guest_path, complex_path,
            model_key=model
        )
        
        if result:
            results.append({
                "model": model,
                "E_interaction_kJ_mol": result['E_interaction_kJ_mol']
            })
        else:
            results.append({
                "model": model,
                "E_interaction_kJ_mol": None
            })
    
    # 汇总
    print(f"\n{'='*60}")
    print("相互作用能比较汇总")
    print(f"{'='*60}")
    print(f"{'模型':<20} {'ΔE (kJ/mol)':<15}")
    print("-" * 35)
    
    for r in results:
        if r['E_interaction_kJ_mol'] is not None:
            print(f"{r['model']:<20} {r['E_interaction_kJ_mol']:<15.2f}")
        else:
            print(f"{r['model']:<20} {'失败'}")
    
    return results


def analyze_adsorption_site_preference(
    mof_path: str,
    c2h2_xyz_content: str,
    sites: dict,
    model_key: str = "mace-mp"
):
    """
    分析不同吸附位点的能量偏好
    
    Args:
        mof_path: MOF 结构文件
        c2h2_xyz_content: 乙炔 XYZ 内容
        sites: 吸附位点字典 {name: (x, y, z)}
        model_key: 模型
    """
    print(f"\n{'='*60}")
    print(f"吸附位点偏好分析")
    print(f"{'='*60}")
    print(f"MOF: {Path(mof_path).stem}")
    print(f"模型: {model_key}")
    print(f"测试位点: {list(sites.keys())}")
    
    # 此功能需要更复杂的实现
    # 这里只展示框架
    print("\n注意: 此功能需要手动准备不同吸附位点的复合物结构")
    print("请参考 structures/HKUST-1_with_C2H2.cif 作为模板")
    
    return None


# ========== 主程序 ==========
if __name__ == "__main__":
    structures_dir = Path(__file__).parent / "structures"
    
    # 示例: 计算 HKUST-1 吸附乙炔的相互作用能
    # 需要三个结构文件：
    # 1. 空的 MOF (宿主)
    # 2. 单独的乙炔分子 (客体)
    # 3. MOF + 乙炔复合物
    
    host_path = structures_dir / "HKUST-1.cif"
    guest_path = structures_dir / "acetylene.xyz"
    complex_path = structures_dir / "HKUST-1_with_C2H2.cif"
    
    if all(p.exists() for p in [host_path, guest_path, complex_path]):
        calculate_interaction_energy_manual(
            str(host_path),
            str(guest_path),
            str(complex_path),
            model_key="mace-mp",
            enable_d3=True
        )
        
        # 多模型比较
        # compare_models_interaction_energy(
        #     str(host_path),
        #     str(guest_path),
        #     str(complex_path),
        #     models=["mace-mp", "orb-v2"]
        # )
    else:
        print("请确保以下文件存在:")
        print(f"  - {host_path}")
        print(f"  - {guest_path}")
        print(f"  - {complex_path}")
