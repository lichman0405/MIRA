# ASE 版本兼容性修复说明

## 问题描述

在 Docker 容器启动时出现错误：
```
ImportError: cannot import name 'ExpCellFilter' from 'ase.constraints'
```

## 根本原因

ASE (Atomic Simulation Environment) 在版本更新过程中进行了重构：

- **ASE < 3.23**: `ExpCellFilter` 和 `StrainFilter` 位于 `ase.constraints` 模块
- **ASE >= 3.23**: 这些类被移到了 `ase.filters` 模块

同时，MIRA 项目需要 **ASE >= 3.27.0** 才能使用 MTKNPT (等温等压系综分子动力学) 等新功能。

## 解决方案

### 1. 修改导入语句（兼容性导入）

在所有使用 `ExpCellFilter` 和 `StrainFilter` 的文件中，使用 try-except 块确保兼容性：

**修改文件：**
- `services/shared/worker_base.py`
- `app/services/optimization.py`

**修改前：**
```python
from ase.constraints import ExpCellFilter, StrainFilter
```

**修改后：**
```python
# ASE 版本兼容性：ExpCellFilter 和 StrainFilter 在新版本中移到了 ase.filters
try:
    from ase.filters import ExpCellFilter, StrainFilter
except ImportError:
    from ase.constraints import ExpCellFilter, StrainFilter
```

这样做的好处：
- ✅ 优先使用新版本 ASE 的导入路径（`ase.filters`）
- ✅ 向后兼容旧版本 ASE（如果有的话）
- ✅ 确保在 ASE >= 3.27.0 环境中正常工作

### 2. 明确 Dockerfile 中的 ASE 版本要求

在所有 Worker Dockerfile 中添加注释说明：

**修改文件：**
- `docker/Dockerfile.mace-orb`
- `docker/Dockerfile.mace-orb-cpu`
- `docker/Dockerfile.matgl`
- `docker/Dockerfile.matgl-cpu`
- （其他 worker Dockerfile 如果安装 ASE）

**修改前：**
```dockerfile
RUN pip install \
    ase>=3.27.0 \
    ...
```

**修改后：**
```dockerfile
# ASE >= 3.27.0 必需：支持 MTKNPT 等新功能，ExpCellFilter/StrainFilter 在 ase.filters 模块
RUN pip install \
    "ase>=3.27.0" \
    ...
```

注意：给 `ase>=3.27.0` 加引号确保版本要求被正确解析。

## 验证测试

创建了测试脚本 `scripts/test_ase_imports.py` 用于验证：

```bash
# 在容器内运行
python scripts/test_ase_imports.py
```

测试内容：
1. ✅ ASE 版本 >= 3.27.0
2. ✅ 成功导入 ExpCellFilter, StrainFilter
3. ✅ NPT 类可用
4. ✅ worker_base.py 正常导入

## 重新构建和部署

### CPU 测试环境

```bash
# 重新构建镜像
./scripts/deploy.sh build-cpu

# 启动服务
./scripts/deploy.sh test-cpu

# 检查日志
docker-compose -f deployments/docker-compose.cpu-test.yml logs mace-orb-worker
```

### GPU 生产环境

```bash
# 重新构建
./scripts/deploy.sh build

# 启动服务
./scripts/deploy.sh up

# 检查所有 worker 日志
docker-compose -f deployments/docker-compose.prod.yml logs
```

## 影响的模块

### 直接影响（需要修改）
- ✅ `services/shared/worker_base.py` - Worker 基类
- ✅ `app/services/optimization.py` - 结构优化服务

### 间接影响（继承自修复的基类）
- `services/mace_orb/` - MACE+ORB Worker
- `services/matgl/` - MatGL Worker
- `services/grace/` - GRACE Worker (如果启用)
- `services/fairchem_sevennet/` - FAIRChem Worker (如果启用)
- `services/mattersim/` - MatterSim Worker (如果启用)

所有这些 Worker 服务都继承自 `BaseWorkerApp`，因此会自动受益于修复。

## ASE 版本历史参考

| ASE 版本 | 发布日期 | 重要变更 |
|---------|---------|---------|
| 3.22.0 | 2022-06 | 稳定版本 |
| 3.23.0 | 2023-11 | `ExpCellFilter` 移到 `ase.filters` |
| 3.27.0 | 2024-04 | 新增 MTKNPT, NPTBerendsen 等 NPT 功能 |
| 3.28.0 | 2024-09 | 当前最新稳定版 |

## 相关资源

- [ASE 官方文档 - Filters](https://wiki.fysik.dtu.dk/ase/ase/filters.html)
- [ASE 官方文档 - NPT](https://wiki.fysik.dtu.dk/ase/ase/md.html#npt-dynamics)
- [ASE GitHub 更新日志](https://gitlab.com/ase/ase/-/blob/master/CHANGELOG.md)

## 总结

这次修复确保了：
1. ✅ MIRA 在 ASE >= 3.27.0 环境下正常工作
2. ✅ 支持所有需要的 NPT 动力学功能
3. ✅ 向后兼容（理论上，虽然我们强制 >= 3.27.0）
4. ✅ 所有 Worker 服务统一使用正确的导入路径

**下次部署时请确保重新构建所有镜像！**
