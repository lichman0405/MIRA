#!/usr/bin/env python3
"""
MIRA - ML 力场模型安装脚本

此脚本帮助安装 MIRA 支持的各种 ML 力场模型。
每个模型都有独立的安装选项，可以按需安装。

使用方法:
    python scripts/install_models.py --all          # 安装所有模型
    python scripts/install_models.py --mace         # 只安装 MACE
    python scripts/install_models.py --mace --orb   # 安装 MACE 和 ORB
    python scripts/install_models.py --check        # 检查已安装的模型
    python scripts/install_models.py --minimal      # 最小安装（仅 MACE）

Author: Shibo Li (shadow.li981@gmail.com)
"""
import subprocess
import sys
import argparse
from typing import List, Dict, Tuple


# ============================================
# 模型安装配置
# ============================================
MODEL_PACKAGES: Dict[str, Dict] = {
    "mace": {
        "name": "MACE",
        "packages": ["mace-torch"],
        "description": "MACE - Multi-ACE models (MP, OFF23, OMAT, MPA, ANI)",
        "models": ["mace-mp", "mace-off23", "mace-omat", "mace-mpa", "mace-ani"]
    },
    "orb": {
        "name": "ORB",
        "packages": ["orb-models"],
        "description": "ORB - Orbital Materials Foundation Model",
        "models": ["orb-v2", "orb-d3-v2", "orb-v3", "orb-v3-mpa"]
    },
    "fairchem": {
        "name": "FAIRChem/OMAT24",
        "packages": ["fairchem-core"],
        "description": "FAIRChem - EquiformerV2, eSEN models",
        "models": ["omat24-base", "omat24-large", "eqv2-omat", "eqv2-mptrj"]
    },
    "grace": {
        "name": "GRACE",
        "packages": ["tensorpotential"],
        "description": "GRACE - Graph-based Atomic Cluster Expansion",
        "models": ["grace-2l", "grace-2l-omat", "grace-2m"]
    },
    "mattersim": {
        "name": "MatterSim",
        "packages": ["mattersim"],
        "description": "MatterSim - Microsoft materials simulator",
        "models": ["mattersim-5m"]
    },
    "sevennet": {
        "name": "SevenNet",
        "packages": ["sevenn"],
        "description": "SevenNet - Neural network interatomic potential",
        "models": ["sevennet-0", "sevennet-mf-ompa", "sevennet-l3i5"]
    },
    "matgl": {
        "name": "MatGL",
        "packages": ["matgl"],
        "description": "MatGL - Materials Graph Library (M3GNet, CHGNet)",
        "models": ["m3gnet", "chgnet"]
    },
}


# ============================================
# 安装函数
# ============================================
def run_pip_install(packages: List[str], upgrade: bool = False) -> bool:
    """运行 pip install"""
    cmd = [sys.executable, "-m", "pip", "install"]
    if upgrade:
        cmd.append("--upgrade")
    cmd.extend(packages)
    
    print(f"  运行: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    return result.returncode == 0


def install_base_dependencies() -> bool:
    """安装基础依赖"""
    print("\n" + "=" * 60)
    print("安装基础依赖...")
    print("=" * 60)
    
    base_packages = [
        "torch>=2.0.0",
        "ase>=3.27.0",
        "phonopy>=2.20.0",
        "torch-dftd>=0.4.0",
        "numpy>=1.24.0",
        "scipy>=1.10.0",
    ]
    
    return run_pip_install(base_packages)


def install_model(model_key: str) -> bool:
    """安装指定模型"""
    if model_key not in MODEL_PACKAGES:
        print(f"  ❌ 未知模型: {model_key}")
        return False
    
    config = MODEL_PACKAGES[model_key]
    print(f"\n📦 安装 {config['name']}...")
    print(f"   {config['description']}")
    
    success = run_pip_install(config["packages"])
    
    if success:
        print(f"  ✅ {config['name']} 安装成功")
    else:
        print(f"  ❌ {config['name']} 安装失败")
    
    return success


def check_model_availability() -> Dict[str, bool]:
    """检查模型可用性"""
    results = {}
    
    print("\n" + "=" * 60)
    print("检查 ML 力场模型安装状态...")
    print("=" * 60 + "\n")
    
    # 检查 MACE
    try:
        from mace.calculators import mace_mp
        results["mace"] = True
        print("✅ MACE        - 已安装")
    except ImportError:
        results["mace"] = False
        print("❌ MACE        - 未安装  (pip install mace-torch)")
    
    # 检查 ORB
    try:
        from orb_models.forcefield import pretrained
        results["orb"] = True
        print("✅ ORB         - 已安装")
    except ImportError:
        results["orb"] = False
        print("❌ ORB         - 未安装  (pip install orb-models)")
    
    # 检查 FAIRChem/OMAT24
    try:
        from fairchem.core import OCPCalculator
        results["fairchem"] = True
        print("✅ FAIRChem    - 已安装")
    except ImportError:
        results["fairchem"] = False
        print("❌ FAIRChem    - 未安装  (pip install fairchem-core)")
    
    # 检查 GRACE
    try:
        import tensorpotential
        results["grace"] = True
        print("✅ GRACE       - 已安装")
    except ImportError:
        results["grace"] = False
        print("❌ GRACE       - 未安装  (pip install tensorpotential)")
    
    # 检查 MatterSim
    try:
        from mattersim.forcefield import MatterSimCalculator
        results["mattersim"] = True
        print("✅ MatterSim   - 已安装")
    except ImportError:
        results["mattersim"] = False
        print("❌ MatterSim   - 未安装  (pip install mattersim)")
    
    # 检查 SevenNet
    try:
        from sevenn.sevennet_calculator import SevenNetCalculator
        results["sevennet"] = True
        print("✅ SevenNet    - 已安装")
    except ImportError:
        results["sevennet"] = False
        print("❌ SevenNet    - 未安装  (pip install sevenn)")
    
    # 检查 MatGL
    try:
        import matgl
        results["matgl"] = True
        print("✅ MatGL       - 已安装")
    except ImportError:
        results["matgl"] = False
        print("❌ MatGL       - 未安装  (pip install matgl)")
    
    # 检查基础依赖
    print("\n--- 基础依赖 ---")
    try:
        import ase
        print(f"✅ ASE         - v{ase.__version__}")
    except ImportError:
        print("❌ ASE         - 未安装")
    
    try:
        import torch
        cuda_status = "CUDA ✓" if torch.cuda.is_available() else "CPU only"
        print(f"✅ PyTorch     - v{torch.__version__} ({cuda_status})")
    except ImportError:
        print("❌ PyTorch     - 未安装")
    
    try:
        import phonopy
        print(f"✅ Phonopy     - v{phonopy.__version__}")
    except ImportError:
        print("❌ Phonopy     - 未安装")
    
    # 统计
    installed = sum(1 for v in results.values() if v)
    total = len(results)
    print(f"\n📊 已安装: {installed}/{total} 个模型家族")
    
    return results


def get_available_models() -> List[str]:
    """获取当前可用的模型列表"""
    available = []
    results = {}
    
    # 静默检查
    try:
        from mace.calculators import mace_mp
        results["mace"] = True
    except ImportError:
        results["mace"] = False
    
    try:
        from orb_models.forcefield import pretrained
        results["orb"] = True
    except ImportError:
        results["orb"] = False
    
    try:
        from fairchem.core import OCPCalculator
        results["fairchem"] = True
    except ImportError:
        results["fairchem"] = False
    
    try:
        import tensorpotential
        results["grace"] = True
    except ImportError:
        results["grace"] = False
    
    try:
        from mattersim.forcefield import MatterSimCalculator
        results["mattersim"] = True
    except ImportError:
        results["mattersim"] = False
    
    try:
        from sevenn.sevennet_calculator import SevenNetCalculator
        results["sevennet"] = True
    except ImportError:
        results["sevennet"] = False
    
    try:
        import matgl
        results["matgl"] = True
    except ImportError:
        results["matgl"] = False
    
    for key, installed in results.items():
        if installed and key in MODEL_PACKAGES:
            available.extend(MODEL_PACKAGES[key]["models"])
    
    return available


# ============================================
# 主函数
# ============================================
def main():
    parser = argparse.ArgumentParser(
        description="MIRA - ML 力场模型安装脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python install_models.py --check              检查已安装的模型
  python install_models.py --minimal            安装最小依赖（仅 MACE）
  python install_models.py --mace --orb         安装 MACE 和 ORB
  python install_models.py --all                安装所有模型
  python install_models.py --recommended        安装推荐模型组合

模型家族:
  --mace       MACE 系列 (mace-mp, mace-off23, mace-omat, ...)
  --orb        ORB 系列 (orb-v2, orb-d3-v2, orb-v3, ...)
  --fairchem   FAIRChem/OMAT24 (omat24-base, eqv2-omat, ...)
  --grace      GRACE 系列 (grace-2l, grace-2m, ...)
  --mattersim  MatterSim (mattersim-5m)
  --sevennet   SevenNet (sevennet-0, sevennet-mf-ompa, ...)
  --matgl      MatGL (m3gnet, chgnet)
"""
    )
    
    # 模型选项
    parser.add_argument("--mace", action="store_true", help="安装 MACE")
    parser.add_argument("--orb", action="store_true", help="安装 ORB")
    parser.add_argument("--fairchem", action="store_true", help="安装 FAIRChem/OMAT24")
    parser.add_argument("--grace", action="store_true", help="安装 GRACE")
    parser.add_argument("--mattersim", action="store_true", help="安装 MatterSim")
    parser.add_argument("--sevennet", action="store_true", help="安装 SevenNet")
    parser.add_argument("--matgl", action="store_true", help="安装 MatGL")
    
    # 快捷选项
    parser.add_argument("--all", action="store_true", help="安装所有模型")
    parser.add_argument("--minimal", action="store_true", help="最小安装（仅 MACE）")
    parser.add_argument("--recommended", action="store_true", 
                       help="推荐安装（MACE + ORB + MatGL）")
    
    # 其他选项
    parser.add_argument("--check", action="store_true", help="检查模型安装状态")
    parser.add_argument("--base", action="store_true", help="只安装基础依赖")
    parser.add_argument("--list", action="store_true", help="列出可用模型")
    
    args = parser.parse_args()
    
    print("""
╔══════════════════════════════════════════════════════════════╗
║                    MIRA Model Installer                       ║
║         MiQroEra Interatomic-potential Reliability Arena      ║
╚══════════════════════════════════════════════════════════════╝
    """)
    
    # 检查模式
    if args.check:
        check_model_availability()
        return 0
    
    # 列出模式
    if args.list:
        available = get_available_models()
        if available:
            print("当前可用的模型:")
            for m in available:
                print(f"  - {m}")
        else:
            print("没有安装任何 ML 力场模型")
            print("运行 `python install_models.py --minimal` 安装基础模型")
        return 0
    
    # 确定要安装的模型
    models_to_install = []
    
    if args.all:
        models_to_install = list(MODEL_PACKAGES.keys())
    elif args.minimal:
        models_to_install = ["mace"]
    elif args.recommended:
        models_to_install = ["mace", "orb", "matgl"]
    else:
        if args.mace:
            models_to_install.append("mace")
        if args.orb:
            models_to_install.append("orb")
        if args.fairchem:
            models_to_install.append("fairchem")
        if args.grace:
            models_to_install.append("grace")
        if args.mattersim:
            models_to_install.append("mattersim")
        if args.sevennet:
            models_to_install.append("sevennet")
        if args.matgl:
            models_to_install.append("matgl")
    
    # 只安装基础依赖
    if args.base:
        install_base_dependencies()
        return 0
    
    # 没有选择任何模型
    if not models_to_install:
        print("请指定要安装的模型，或使用 --help 查看帮助")
        print("\n快速开始:")
        print("  python install_models.py --minimal      # 安装 MACE")
        print("  python install_models.py --recommended  # 推荐组合")
        print("  python install_models.py --check        # 检查状态")
        return 1
    
    # 安装基础依赖
    install_base_dependencies()
    
    # 安装选定的模型
    print("\n" + "=" * 60)
    print(f"将安装以下模型: {', '.join(models_to_install)}")
    print("=" * 60)
    
    success_count = 0
    failed = []
    
    for model in models_to_install:
        if install_model(model):
            success_count += 1
        else:
            failed.append(model)
    
    # 总结
    print("\n" + "=" * 60)
    print("安装完成!")
    print("=" * 60)
    print(f"✅ 成功: {success_count}/{len(models_to_install)}")
    
    if failed:
        print(f"❌ 失败: {', '.join(failed)}")
        print("\n失败的模型可能需要手动安装，请参考各模型的官方文档")
    
    # 最终检查
    print("\n运行检查...")
    check_model_availability()
    
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
