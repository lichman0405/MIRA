#!/usr/bin/env python3
"""
测试 ASE 导入兼容性
验证 ExpCellFilter 和 StrainFilter 的导入位置
"""
import sys

def test_ase_version():
    """测试 ASE 版本"""
    try:
        import ase
        print(f"✓ ASE 版本: {ase.__version__}")
        
        # 检查版本是否满足要求
        from packaging import version
        ase_version = version.parse(ase.__version__)
        required_version = version.parse("3.27.0")
        
        if ase_version >= required_version:
            print(f"✓ ASE 版本满足要求 (>= 3.27.0)")
        else:
            print(f"✗ ASE 版本过低: {ase.__version__} < 3.27.0")
            return False
            
    except ImportError as e:
        print(f"✗ ASE 未安装: {e}")
        return False
    
    return True


def test_filters_import():
    """测试 filters 导入"""
    print("\n=== 测试 Filters 导入 ===")
    
    # 测试新版本导入路径 (ase.filters)
    try:
        from ase.filters import ExpCellFilter, StrainFilter
        print("✓ 从 ase.filters 成功导入 ExpCellFilter, StrainFilter")
        print(f"  ExpCellFilter: {ExpCellFilter}")
        print(f"  StrainFilter: {StrainFilter}")
        return True
    except ImportError as e:
        print(f"✗ 从 ase.filters 导入失败: {e}")
    
    # 测试旧版本导入路径 (ase.constraints)
    try:
        from ase.constraints import ExpCellFilter, StrainFilter
        print("✓ 从 ase.constraints 成功导入 ExpCellFilter, StrainFilter (旧版本)")
        print(f"  ExpCellFilter: {ExpCellFilter}")
        print(f"  StrainFilter: {StrainFilter}")
        return True
    except ImportError as e:
        print(f"✗ 从 ase.constraints 导入失败: {e}")
    
    return False


def test_npt_classes():
    """测试 NPT 相关类"""
    print("\n=== 测试 NPT 类 ===")
    
    npt_classes = [
        ("ase.md.npt", "NPT"),
        ("ase.md.nptberendsen", "NPTBerendsen"),
        ("ase.md.npt", "NPTBerendsen"),  # 新版本可能在这里
    ]
    
    found_npt = False
    for module_name, class_name in npt_classes:
        try:
            module = __import__(module_name, fromlist=[class_name])
            cls = getattr(module, class_name)
            print(f"✓ 找到 {module_name}.{class_name}")
            found_npt = True
        except (ImportError, AttributeError):
            pass
    
    if not found_npt:
        print("✗ 未找到 NPT 相关类")
        return False
    
    return True


def test_worker_base_import():
    """测试 worker_base.py 导入"""
    print("\n=== 测试 worker_base.py 导入 ===")
    
    try:
        # 修改 sys.path 以便导入
        import os
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'services'))
        
        from shared.worker_base import BaseWorkerApp
        print("✓ 成功导入 BaseWorkerApp")
        return True
    except ImportError as e:
        print(f"✗ 导入 BaseWorkerApp 失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """主函数"""
    print("=" * 60)
    print("ASE 兼容性测试")
    print("=" * 60)
    
    results = []
    
    # 测试 ASE 版本
    results.append(("ASE 版本", test_ase_version()))
    
    # 测试 Filters 导入
    results.append(("Filters 导入", test_filters_import()))
    
    # 测试 NPT 类
    results.append(("NPT 类", test_npt_classes()))
    
    # 测试 worker_base 导入
    results.append(("worker_base 导入", test_worker_base_import()))
    
    # 总结
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)
    
    all_passed = True
    for name, passed in results:
        status = "✓ 通过" if passed else "✗ 失败"
        print(f"{name:20s} {status}")
        if not passed:
            all_passed = False
    
    print("=" * 60)
    
    if all_passed:
        print("✓ 所有测试通过！")
        return 0
    else:
        print("✗ 部分测试失败")
        return 1


if __name__ == "__main__":
    sys.exit(main())
