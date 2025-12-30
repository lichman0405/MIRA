#!/usr/bin/env python3
"""
修复 PyTorch 和 torchvision 版本兼容性问题

错误: RuntimeError: operator torchvision::nms does not exist

原因: PyTorch 和 torchvision 版本不匹配
"""
import subprocess
import sys


def get_current_versions():
    """获取当前安装的包版本"""
    print("=" * 60)
    print("检查当前版本")
    print("=" * 60)
    
    packages = ["torch", "torchvision", "torchaudio", "torchmetrics", "mace-torch"]
    
    for pkg in packages:
        try:
            result = subprocess.run(
                [sys.executable, "-m", "pip", "show", pkg],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                for line in result.stdout.split("\n"):
                    if line.startswith("Version:"):
                        version = line.split(":")[1].strip()
                        print(f"  {pkg:20s} {version}")
                        break
            else:
                print(f"  {pkg:20s} 未安装")
        except Exception as e:
            print(f"  {pkg:20s} 检查失败: {e}")
    
    print()


def uninstall_packages():
    """卸载有问题的包"""
    print("=" * 60)
    print("卸载旧版本包")
    print("=" * 60)
    
    # 卸载顺序很重要：先卸载依赖这些包的，再卸载基础包
    packages_to_uninstall = [
        "mace-torch",
        "torchmetrics",
        "torchvision",
        "torchaudio",
        "torch"
    ]
    
    for pkg in packages_to_uninstall:
        print(f"\n卸载 {pkg}...")
        subprocess.run(
            [sys.executable, "-m", "pip", "uninstall", "-y", pkg],
            capture_output=True
        )
    
    print("\n✓ 卸载完成")


def install_compatible_versions():
    """安装兼容的版本"""
    print("\n" + "=" * 60)
    print("安装兼容版本")
    print("=" * 60)
    
    # PyTorch 2.5.1 + torchvision 0.20.1 (CUDA 12.4)
    # 或 PyTorch 2.4.1 + torchvision 0.19.1 (CUDA 12.1)
    
    print("\n选择安装版本:")
    print("1. PyTorch 2.5.1 + CUDA 12.4 (推荐，最新)")
    print("2. PyTorch 2.4.1 + CUDA 12.1 (稳定)")
    print("3. CPU 版本 (无 GPU)")
    
    choice = input("\n请选择 (1/2/3): ").strip()
    
    if choice == "1":
        print("\n安装 PyTorch 2.5.1 + CUDA 12.4...")
        cmd = [
            sys.executable, "-m", "pip", "install",
            "torch==2.5.1",
            "torchvision==0.20.1",
            "torchaudio==2.5.1",
            "--index-url", "https://download.pytorch.org/whl/cu124"
        ]
    elif choice == "2":
        print("\n安装 PyTorch 2.4.1 + CUDA 12.1...")
        cmd = [
            sys.executable, "-m", "pip", "install",
            "torch==2.4.1",
            "torchvision==0.19.1",
            "torchaudio==2.4.1",
            "--index-url", "https://download.pytorch.org/whl/cu121"
        ]
    elif choice == "3":
        print("\n安装 PyTorch CPU 版本...")
        cmd = [
            sys.executable, "-m", "pip", "install",
            "torch",
            "torchvision",
            "torchaudio",
            "--index-url", "https://download.pytorch.org/whl/cpu"
        ]
    else:
        print("无效选择，退出")
        return False
    
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("✗ PyTorch 安装失败")
        return False
    
    print("\n✓ PyTorch 安装成功")
    
    # 安装兼容的 torchmetrics
    print("\n安装 torchmetrics...")
    subprocess.run([
        sys.executable, "-m", "pip", "install",
        "torchmetrics>=1.0.0"
    ])
    
    # 重新安装 mace-torch
    print("\n重新安装 mace-torch...")
    subprocess.run([
        sys.executable, "-m", "pip", "install",
        "mace-torch"
    ])
    
    print("\n✓ 所有包安装完成")
    return True


def verify_installation():
    """验证安装"""
    print("\n" + "=" * 60)
    print("验证安装")
    print("=" * 60)
    
    get_current_versions()
    
    print("测试导入...")
    test_imports = [
        "import torch",
        "import torchvision",
        "from mace.calculators import mace_mp"
    ]
    
    all_passed = True
    for import_statement in test_imports:
        try:
            exec(import_statement)
            module_name = import_statement.split()[1].split(".")[0]
            print(f"  ✓ {module_name:20s} 导入成功")
        except Exception as e:
            module_name = import_statement.split()[1].split(".")[0]
            print(f"  ✗ {module_name:20s} 导入失败: {e}")
            all_passed = False
    
    return all_passed


def main():
    """主函数"""
    print("\n" + "=" * 60)
    print("PyTorch/torchvision 兼容性修复工具")
    print("=" * 60)
    
    # 显示当前版本
    get_current_versions()
    
    # 询问是否继续
    response = input("\n是否继续修复? (y/n): ").strip().lower()
    if response != 'y':
        print("取消操作")
        return
    
    # 卸载旧版本
    uninstall_packages()
    
    # 安装兼容版本
    if not install_compatible_versions():
        print("\n安装失败，请手动修复")
        return
    
    # 验证
    if verify_installation():
        print("\n" + "=" * 60)
        print("✓ 修复成功！")
        print("=" * 60)
        print("\n现在可以运行 examples 了:")
        print("  python examples/01_basic_usage.py")
    else:
        print("\n" + "=" * 60)
        print("⚠ 部分导入仍有问题，请检查错误信息")
        print("=" * 60)


if __name__ == "__main__":
    main()
