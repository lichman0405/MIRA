"""
MIRA Examples - 依赖检查与模型安装帮助模块

在运行 examples 之前，先检查必要的依赖是否已安装。
如果缺少依赖，会提供安装指导。

使用方法:
    from setup_check import ensure_dependencies, check_available_models
    
    # 检查基础依赖
    ensure_dependencies()
    
    # 获取可用模型
    models = check_available_models()
"""
import sys
import subprocess
from typing import List, Dict, Tuple
from pathlib import Path


# ============================================
# ANSI 颜色代码
# ============================================
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

    @classmethod
    def disable(cls):
        """禁用颜色（Windows 兼容）"""
        cls.GREEN = cls.RED = cls.YELLOW = cls.BLUE = cls.BOLD = cls.END = ''


# Windows 兼容
if sys.platform == 'win32':
    try:
        import colorama
        colorama.init()
    except ImportError:
        Colors.disable()


# ============================================
# 模型配置
# ============================================
MODEL_REGISTRY = {
    # MACE 系列
    "mace-mp": ("mace", "from mace.calculators import mace_mp"),
    "mace-off23": ("mace", "from mace.calculators import mace_mp"),
    "mace-omat": ("mace", "from mace.calculators import MACECalculator"),
    "mace-mpa": ("mace", "from mace.calculators import mace_mp"),
    "mace-ani": ("mace", "from mace.calculators import mace_mp"),
    "mace_prod_b3": ("mace", "from mace.calculators import mace_mp"),
    "mace_prod": ("mace", "from mace.calculators import mace_mp"),
    
    # ORB 系列
    "orb-v2": ("orb", "from orb_models.forcefield import pretrained"),
    "orb-d3-v2": ("orb", "from orb_models.forcefield import pretrained"),
    "orb-v3": ("orb", "from orb_models.forcefield import pretrained"),
    "orb-v3-mpa": ("orb", "from orb_models.forcefield import pretrained"),
    "orb_prod": ("orb", "from orb_models.forcefield import pretrained"),
    
    # FAIRChem/OMAT24
    "omat24-base": ("fairchem", "from fairchem.core import OCPCalculator"),
    "omat24-large": ("fairchem", "from fairchem.core import OCPCalculator"),
    "eqv2-omat": ("fairchem", "from fairchem.core import OCPCalculator"),
    "eqv2-mptrj": ("fairchem", "from fairchem.core import OCPCalculator"),
    
    # GRACE
    "grace-2l": ("grace", "import tensorpotential"),
    "grace-2l-omat": ("grace", "import tensorpotential"),
    "grace-2m": ("grace", "import tensorpotential"),
    
    # MatterSim
    "mattersim-5m": ("mattersim", "from mattersim.forcefield import MatterSimCalculator"),
    
    # SevenNet
    "sevennet-0": ("sevennet", "from sevenn.sevennet_calculator import SevenNetCalculator"),
    "sevennet-mf-ompa": ("sevennet", "from sevenn.sevennet_calculator import SevenNetCalculator"),
    "sevennet-l3i5": ("sevennet", "from sevenn.sevennet_calculator import SevenNetCalculator"),
    
    # MatGL
    "m3gnet": ("matgl", "import matgl"),
    "chgnet": ("matgl", "import matgl"),
}

FAMILY_INSTALL_CMD = {
    "mace": "pip install mace-torch",
    "orb": "pip install orb-models",
    "fairchem": "pip install fairchem-core",
    "grace": "pip install tensorpotential",
    "mattersim": "pip install mattersim",
    "sevennet": "pip install sevenn",
    "matgl": "pip install matgl",
}


# ============================================
# 检查函数
# ============================================
def check_package(import_statement: str) -> bool:
    """检查包是否可导入"""
    try:
        exec(import_statement)
        return True
    except ImportError:
        return False


def check_base_dependencies() -> Dict[str, Tuple[bool, str]]:
    """检查基础依赖"""
    results = {}
    
    # ASE
    try:
        import ase
        results["ase"] = (True, ase.__version__)
    except ImportError:
        results["ase"] = (False, "未安装 - pip install ase>=3.27.0")
    
    # PyTorch
    try:
        import torch
        cuda_info = "CUDA ✓" if torch.cuda.is_available() else "CPU"
        results["torch"] = (True, f"{torch.__version__} ({cuda_info})")
    except ImportError:
        results["torch"] = (False, "未安装 - pip install torch")
    
    # Phonopy
    try:
        import phonopy
        results["phonopy"] = (True, phonopy.__version__)
    except ImportError:
        results["phonopy"] = (False, "未安装 - pip install phonopy")
    
    # NumPy
    try:
        import numpy
        results["numpy"] = (True, numpy.__version__)
    except ImportError:
        results["numpy"] = (False, "未安装 - pip install numpy")
    
    # Requests
    try:
        import requests
        results["requests"] = (True, requests.__version__)
    except ImportError:
        results["requests"] = (False, "未安装 - pip install requests")
    
    return results


def check_model_families() -> Dict[str, bool]:
    """检查各模型家族的可用性"""
    results = {}
    
    # MACE
    results["mace"] = check_package("from mace.calculators import mace_mp")
    
    # ORB
    results["orb"] = check_package("from orb_models.forcefield import pretrained")
    
    # FAIRChem
    results["fairchem"] = check_package("from fairchem.core import OCPCalculator")
    
    # GRACE
    results["grace"] = check_package("import tensorpotential")
    
    # MatterSim
    results["mattersim"] = check_package("from mattersim.forcefield import MatterSimCalculator")
    
    # SevenNet
    results["sevennet"] = check_package("from sevenn.sevennet_calculator import SevenNetCalculator")
    
    # MatGL
    results["matgl"] = check_package("import matgl")
    
    return results


def get_available_models() -> List[str]:
    """获取当前可用的模型列表"""
    families = check_model_families()
    available = []
    
    for model_key, (family, _) in MODEL_REGISTRY.items():
        if families.get(family, False):
            available.append(model_key)
    
    return available


def is_model_available(model_key: str) -> bool:
    """检查特定模型是否可用"""
    if model_key not in MODEL_REGISTRY:
        return False
    family, _ = MODEL_REGISTRY[model_key]
    families = check_model_families()
    return families.get(family, False)


# ============================================
# 主要接口
# ============================================
def ensure_dependencies(required_models: List[str] = None, verbose: bool = True):
    """
    确保依赖已安装
    
    Args:
        required_models: 需要的模型列表（可选）
        verbose: 是否打印详细信息
    """
    if verbose:
        print(f"\n{Colors.BOLD}{'='*60}{Colors.END}")
        print(f"{Colors.BOLD}MIRA 依赖检查{Colors.END}")
        print(f"{'='*60}\n")
    
    # 检查基础依赖
    base_deps = check_base_dependencies()
    missing_base = []
    
    if verbose:
        print(f"{Colors.BOLD}基础依赖:{Colors.END}")
    
    for name, (installed, info) in base_deps.items():
        if installed:
            if verbose:
                print(f"  {Colors.GREEN}✓{Colors.END} {name}: {info}")
        else:
            missing_base.append(name)
            if verbose:
                print(f"  {Colors.RED}✗{Colors.END} {name}: {info}")
    
    if missing_base:
        print(f"\n{Colors.RED}错误: 缺少基础依赖{Colors.END}")
        print("请运行: pip install ase>=3.27.0 torch phonopy numpy requests")
        sys.exit(1)
    
    # 检查 ML 模型
    if verbose:
        print(f"\n{Colors.BOLD}ML 力场模型:{Colors.END}")
    
    families = check_model_families()
    installed_families = []
    
    for family, installed in families.items():
        if installed:
            installed_families.append(family)
            if verbose:
                print(f"  {Colors.GREEN}✓{Colors.END} {family.upper()}")
        else:
            if verbose:
                print(f"  {Colors.YELLOW}○{Colors.END} {family.upper()} (未安装)")
    
    # 检查特定模型要求
    if required_models:
        missing_models = []
        for model in required_models:
            if not is_model_available(model):
                missing_models.append(model)
        
        if missing_models:
            print(f"\n{Colors.YELLOW}警告: 以下模型不可用: {', '.join(missing_models)}{Colors.END}")
            
            # 找出需要安装的包
            needed_families = set()
            for model in missing_models:
                if model in MODEL_REGISTRY:
                    family, _ = MODEL_REGISTRY[model]
                    needed_families.add(family)
            
            print("\n请安装以下包:")
            for family in needed_families:
                print(f"  {FAMILY_INSTALL_CMD.get(family, f'pip install {family}')}")
    
    if verbose:
        available = get_available_models()
        print(f"\n{Colors.BOLD}可用模型数量: {len(available)}{Colors.END}")
        
        if not installed_families:
            print(f"\n{Colors.YELLOW}提示: 没有安装任何 ML 力场模型{Colors.END}")
            print("运行安装脚本: python scripts/install_models.py --minimal")
    
    if verbose:
        print(f"\n{'='*60}\n")
    
    return len(installed_families) > 0


def check_available_models(verbose: bool = False) -> List[str]:
    """
    获取可用模型列表
    
    Args:
        verbose: 是否打印详细信息
        
    Returns:
        可用模型的键列表
    """
    available = get_available_models()
    
    if verbose:
        print(f"\n{Colors.BOLD}当前可用的模型 ({len(available)} 个):{Colors.END}")
        for model in available:
            print(f"  - {model}")
        print()
    
    return available


def get_first_available_model(preferred: List[str] = None) -> str:
    """
    获取第一个可用的模型
    
    Args:
        preferred: 优先选择的模型列表
        
    Returns:
        模型键，如果没有可用模型返回 None
    """
    available = get_available_models()
    
    if not available:
        return None
    
    if preferred:
        for model in preferred:
            if model in available:
                return model
    
    return available[0]


def print_install_guide():
    """打印安装指南"""
    print(f"""
{Colors.BOLD}{'='*60}
MIRA - ML 力场模型安装指南
{'='*60}{Colors.END}

运行安装脚本:
  python scripts/install_models.py --check        # 检查状态
  python scripts/install_models.py --minimal      # 最小安装 (MACE)
  python scripts/install_models.py --recommended  # 推荐 (MACE+ORB+MatGL)
  python scripts/install_models.py --all          # 全部安装

或手动安装:
  pip install mace-torch      # MACE
  pip install orb-models      # ORB  
  pip install fairchem-core   # FAIRChem/OMAT24
  pip install tensorpotential # GRACE
  pip install mattersim       # MatterSim
  pip install sevenn          # SevenNet
  pip install matgl           # MatGL (M3GNet, CHGNet)

{'='*60}
""")


# ============================================
# 命令行入口
# ============================================
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="MIRA 依赖检查")
    parser.add_argument("--guide", action="store_true", help="显示安装指南")
    parser.add_argument("--list", action="store_true", help="列出可用模型")
    args = parser.parse_args()
    
    if args.guide:
        print_install_guide()
    elif args.list:
        check_available_models(verbose=True)
    else:
        ensure_dependencies(verbose=True)
