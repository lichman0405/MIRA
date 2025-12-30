#!/usr/bin/env python3
"""
MIRA - ML 力场模型安装脚本

██╗    ██╗ █████╗ ██████╗ ███╗   ██╗██╗███╗   ██╗ ██████╗ 
██║    ██║██╔══██╗██╔══██╗████╗  ██║██║████╗  ██║██╔════╝ 
██║ █╗ ██║███████║██████╔╝██╔██╗ ██║██║██╔██╗ ██║██║  ███╗
██║███╗██║██╔══██║██╔══██╗██║╚██╗██║██║██║╚██╗██║██║   ██║
╚███╔███╔╝██║  ██║██║  ██║██║ ╚████║██║██║ ╚████║╚██████╔╝
 ╚══╝╚══╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝╚═╝  ╚═══╝ ╚═════╝ 

⚠️  重要提示：不同 ML 力场模型有不兼容的依赖版本！
    不同模型对 PyTorch、e3nn、TensorFlow 等有不同版本要求，
    在同一环境安装所有模型会导致依赖冲突！

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
推荐策略：使用多 conda 环境，每个环境安装一组兼容的模型
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

兼容的模型组合：

  🅰 组合 A (推荐入门):  MACE + ORB
     依赖: PyTorch + e3nn==0.4.4
     命令: python scripts/install_models.py --combo-a
     适合: 初学者、MOF 基准测试
     
  🅱 组合 B:  FAIRChem (OMAT24) + SevenNet  
     依赖: PyTorch + e3nn>=0.5.0
     命令: python scripts/install_models.py --combo-b
     适合: 大规模材料预测
     
  🅲 组合 C:  MatGL (M3GNet + CHGNet)
     依赖: PyTorch + DGL
     命令: python scripts/install_models.py --combo-c
     适合: 电池材料、晶体结构
     
  🅳 组合 D:  GRACE
     依赖: TensorFlow (与 PyTorch 模型隔离)
     命令: python scripts/install_models.py --combo-d
     适合: 高精度力场

多环境设置示例：
  # 环境 1: MACE + ORB
  conda create -n mira-mace python=3.10 && conda activate mira-mace
  pip install mace-torch orb-models ase phonopy

  # 环境 2: FAIRChem + SevenNet
  conda create -n mira-fairchem python=3.10 && conda activate mira-fairchem  
  pip install fairchem-core sevenn ase phonopy

使用方法:
    python scripts/install_models.py --check           # 检查已安装的模型
    python scripts/install_models.py --combo-a         # 安装组合A (MACE+ORB)
    python scripts/install_models.py --mace            # 只安装 MACE
    python scripts/install_models.py --mace --orb      # 安装 MACE 和 ORB

⚠️ 不推荐使用 --all，会导致依赖冲突！

Author: Shibo Li (shadow.li981@gmail.com)
"""
import subprocess
import sys
import argparse
from typing import List, Dict, Tuple


# ============================================
# 兼容性组合配置
# ============================================
COMPATIBLE_COMBOS = {
    "combo-a": {
        "name": "组合 A (MACE + ORB)",
        "models": ["mace", "orb"],
        "description": "最稳定的组合，适合入门和 MOF 基准测试",
        "deps": "PyTorch + e3nn==0.4.4"
    },
    "combo-b": {
        "name": "组合 B (FAIRChem + SevenNet)",
        "models": ["fairchem", "sevennet"],
        "description": "OMAT24 生态，适合大规模材料预测",
        "deps": "PyTorch + e3nn>=0.5.0"
    },
    "combo-c": {
        "name": "组合 C (MatGL)",
        "models": ["matgl"],
        "description": "M3GNet + CHGNet，适合电池材料、晶体结构",
        "deps": "PyTorch + DGL"
    },
    "combo-d": {
        "name": "组合 D (GRACE)",
        "models": ["grace"],
        "description": "Graph-based ACE，需要 TensorFlow",
        "deps": "TensorFlow"
    },
    "combo-e": {
        "name": "组合 E (MatterSim)",
        "models": ["mattersim"],
        "description": "Microsoft MatterSim",
        "deps": "PyTorch"
    },
}


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
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
⚠️  重要：不同模型有不兼容的依赖，请使用兼容组合安装！
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

兼容组合 (推荐):
  --combo-a     🅰 MACE + ORB           # 推荐入门，最稳定
  --combo-b     🅱 FAIRChem + SevenNet  # OMAT24 生态
  --combo-c     🅲 MatGL (M3GNet+CHGNet)# 独立环境
  --combo-d     🅳 GRACE                # TensorFlow
  --combo-e     🅴 MatterSim            # Microsoft

单独安装:
  --mace        MACE 系列 (mace-mp, mace-off23, mace-omat, ...)
  --orb         ORB 系列 (orb-v2, orb-d3-v2, orb-v3, ...)
  --fairchem    FAIRChem/OMAT24 (omat24-base, eqv2-omat, ...)
  --grace       GRACE 系列 (grace-2l, grace-2m, ...)
  --mattersim   MatterSim (mattersim-5m)
  --sevennet    SevenNet (sevennet-0, sevennet-mf-ompa, ...)
  --matgl       MatGL (m3gnet, chgnet)

示例:
  python install_models.py --check              # 检查已安装的模型
  python install_models.py --combo-a            # 安装组合A (推荐)
  python install_models.py --mace               # 只安装 MACE
  python install_models.py --mace --orb         # 安装 MACE 和 ORB

多环境策略 (推荐):
  conda create -n mira-mace python=3.10
  conda activate mira-mace
  python install_models.py --combo-a
"""
    )
    
    # 兼容组合选项 (推荐)
    combo_group = parser.add_argument_group("兼容组合 (推荐)")
    combo_group.add_argument("--combo-a", action="store_true", 
                            help="🅰 安装 MACE + ORB (推荐入门)")
    combo_group.add_argument("--combo-b", action="store_true", 
                            help="🅱 安装 FAIRChem + SevenNet")
    combo_group.add_argument("--combo-c", action="store_true", 
                            help="🅲 安装 MatGL (M3GNet + CHGNet)")
    combo_group.add_argument("--combo-d", action="store_true", 
                            help="🅳 安装 GRACE (需要 TensorFlow)")
    combo_group.add_argument("--combo-e", action="store_true", 
                            help="🅴 安装 MatterSim")
    
    # 单独模型选项
    model_group = parser.add_argument_group("单独模型")
    model_group.add_argument("--mace", action="store_true", help="安装 MACE")
    model_group.add_argument("--orb", action="store_true", help="安装 ORB")
    model_group.add_argument("--fairchem", action="store_true", help="安装 FAIRChem/OMAT24")
    model_group.add_argument("--grace", action="store_true", help="安装 GRACE")
    model_group.add_argument("--mattersim", action="store_true", help="安装 MatterSim")
    model_group.add_argument("--sevennet", action="store_true", help="安装 SevenNet")
    model_group.add_argument("--matgl", action="store_true", help="安装 MatGL")
    
    # 快捷选项 (保留但添加警告)
    shortcut_group = parser.add_argument_group("快捷选项")
    shortcut_group.add_argument("--all", action="store_true", 
                               help="⚠️ 安装所有模型 (不推荐，会有依赖冲突)")
    shortcut_group.add_argument("--minimal", action="store_true", 
                               help="最小安装（仅 MACE）")
    shortcut_group.add_argument("--recommended", action="store_true", 
                               help="同 --combo-a")
    
    # 其他选项
    other_group = parser.add_argument_group("其他选项")
    other_group.add_argument("--check", action="store_true", help="检查模型安装状态")
    other_group.add_argument("--base", action="store_true", help="只安装基础依赖")
    other_group.add_argument("--list", action="store_true", help="列出可用模型")
    other_group.add_argument("--force", action="store_true", 
                            help="强制安装 (忽略依赖冲突警告)")
    
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
            print("运行 `python install_models.py --combo-a` 安装推荐组合")
        return 0
    
    # 确定要安装的模型
    models_to_install = []
    combo_name = None
    
    # 处理兼容组合
    if args.combo_a or args.recommended:
        combo = COMPATIBLE_COMBOS["combo-a"]
        models_to_install = combo["models"]
        combo_name = combo["name"]
    elif args.combo_b:
        combo = COMPATIBLE_COMBOS["combo-b"]
        models_to_install = combo["models"]
        combo_name = combo["name"]
    elif args.combo_c:
        combo = COMPATIBLE_COMBOS["combo-c"]
        models_to_install = combo["models"]
        combo_name = combo["name"]
    elif args.combo_d:
        combo = COMPATIBLE_COMBOS["combo-d"]
        models_to_install = combo["models"]
        combo_name = combo["name"]
    elif args.combo_e:
        combo = COMPATIBLE_COMBOS["combo-e"]
        models_to_install = combo["models"]
        combo_name = combo["name"]
    elif args.all:
        # --all 需要警告
        if not args.force:
            print("""
╔══════════════════════════════════════════════════════════════╗
║  ⚠️  警告：--all 选项会导致依赖冲突！                         ║
╠══════════════════════════════════════════════════════════════╣
║                                                              ║
║  不同 ML 力场模型有不兼容的依赖版本：                         ║
║  • mace-torch 需要 e3nn==0.4.4                               ║
║  • fairchem-core 需要 e3nn>=0.5.0                            ║
║  • matgl 需要特定版本的 PyTorch                              ║
║  • GRACE 需要 TensorFlow                                     ║
║                                                              ║
║  推荐方案：使用多 conda 环境，每个环境安装一组兼容的模型       ║
║                                                              ║
║  兼容组合：                                                   ║
║    --combo-a    MACE + ORB (推荐入门)                        ║
║    --combo-b    FAIRChem + SevenNet                          ║
║    --combo-c    MatGL (M3GNet + CHGNet)                      ║
║    --combo-d    GRACE                                        ║
║                                                              ║
║  如果仍要继续，请添加 --force 参数：                          ║
║    python install_models.py --all --force                    ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
            """)
            return 1
        else:
            print("""
⚠️  警告：正在安装所有模型，可能会有依赖冲突！
    某些模型可能无法正常工作。
            """)
            models_to_install = list(MODEL_PACKAGES.keys())
    elif args.minimal:
        models_to_install = ["mace"]
    else:
        # 单独模型选项
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
        print("""
请指定要安装的模型组合：

推荐组合 (兼容):
  python install_models.py --combo-a     # 🅰 MACE + ORB (推荐入门)
  python install_models.py --combo-b     # 🅱 FAIRChem + SevenNet
  python install_models.py --combo-c     # 🅲 MatGL
  python install_models.py --combo-d     # 🅳 GRACE
  python install_models.py --combo-e     # 🅴 MatterSim

其他选项:
  python install_models.py --check       # 检查已安装的模型
  python install_models.py --help        # 查看完整帮助
""")
        return 1
    
    # 安装基础依赖
    install_base_dependencies()
    
    # 安装选定的模型
    print("\n" + "=" * 60)
    if combo_name:
        print(f"安装 {combo_name}")
    print(f"包含模型: {', '.join(models_to_install)}")
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
