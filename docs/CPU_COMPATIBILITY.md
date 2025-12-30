# MIRA Examples - CPU 模式兼容性说明

## 概述

所有 MIRA examples 都是通过 HTTP API 调用服务的，因此**理论上可以在任何环境下运行**（本地或远程）。但是，**实际能否成功执行取决于后端服务支持的模型**。

## CPU 模式支持的模型

| 模型家族 | 模型 | CPU 支持 | 说明 |
|----------|------|----------|------|
| **MACE** | mace-mp | ✅ | 完全支持 |
| | mace-off23 | ✅ | 完全支持 |
| | mace-omat | ✅ | 完全支持 |
| | mace-mpa | ✅ | 支持，但速度很慢 |
| | mace-ani | ✅ | 完全支持 |
| **ORB** | orb-v2 | ✅ | 完全支持 |
| | orb-d3-v2 | ✅ | 完全支持 |
| | orb-v3 | ✅ | 完全支持 |
| | orb-v3-mpa | ✅ | 支持，但速度很慢 |
| | orb-omat-v3-lora | ✅ | 完全支持 |
| **MatGL** | m3gnet | ✅ | 完全支持 |
| | chgnet | ✅ | 完全支持 |
| **OMAT24** | omat24-* | ❌ | 仅 GPU，需 FAIRChem Worker |
| | eqv2-* | ❌ | 仅 GPU，需 FAIRChem Worker |
| **SevenNet** | sevennet-* | ❌ | 仅 GPU，需 SevenNet Worker |
| **GRACE** | grace-* | ❌ | 仅 GPU，需 GRACE Worker |
| **MatterSim** | mattersim-* | ❌ | 仅 GPU，需 MatterSim Worker |
| **PosEGNN** | posegnn | ❌ | 仅 GPU，需 FAIRChem Worker |

## Examples 兼容性分析

### ✅ 完全兼容 (CPU 模式下可运行)

以下 examples 默认使用 CPU 兼容的模型，可以在 CPU-only 环境下直接运行：

#### 1. `01_basic_usage.py` - 基础使用
- ✅ **完全兼容**
- 主要测试 API 连接和基本操作
- 不涉及具体模型计算

#### 2. `03_md_stability.py` - MD 稳定性测试
- ✅ **默认兼容**
- 默认使用 `mace-mp`
- CPU 模式下可运行，但速度较慢（MD 步数多）

#### 3. `04_bulk_modulus.py` - 体积模量计算
- ✅ **默认兼容**
- 默认使用 `mace-mp`
- CPU 模式下可运行

#### 4. `05_heat_capacity.py` - 热容计算
- ✅ **默认兼容**
- 默认使用 `mace-mp`
- CPU 模式下可运行，但 Phonopy 计算较慢

#### 5. `06_acetylene_adsorption.py` - 乙炔吸附
- ✅ **默认兼容**
- 默认使用 `mace-mp`
- CPU 模式下可运行

---

### ⚠️ 部分兼容 (需修改配置)

以下 examples 默认测试多个模型，在 CPU 模式下需要修改配置：

#### 6. `02_structure_optimization.py` - 结构优化
- ⚠️ **需要修改**
- 默认测试模型列表：
  ```python
  models = [
      "mace-mp",      # ✅ CPU 支持
      "orb-v2",       # ✅ CPU 支持
      "omat24-base",  # ❌ 仅 GPU
      "grace-2l",     # ❌ 仅 GPU
      "mattersim-5m", # ❌ 仅 GPU
      "sevennet-0"    # ❌ 仅 GPU
  ]
  ```

**CPU 模式修改方法：**
```python
# 修改前（第 164 行左右）
models = [
    "mace-mp",
    "orb-v2", 
    "omat24-base",   # 删除此行
    "grace-2l",      # 删除此行
    "mattersim-5m",  # 删除此行
    "sevennet-0"     # 删除此行
]

# 修改后
models = [
    "mace-mp",
    "orb-v2",
    "m3gnet",       # 添加 MatGL 模型
    "chgnet"        # 添加 CHGNet 模型
]
```

或运行时通过环境变量指定：
```bash
# 只测试 MACE 和 ORB
export MIRA_TEST_MODELS="mace-mp,orb-v2,m3gnet"
python examples/02_structure_optimization.py
```

---

#### 7. `07_full_benchmark.py` - 完整基准测试
- ⚠️ **需要修改**
- 默认测试 20+ 个模型，包含大量 GPU-only 模型

**CPU 模式修改方法：**
```python
# 修改前（第 36-46 行左右）
ALL_MODELS = [
    # MACE 系列
    "mace-mp", "mace-off23", "mace-omat", "mace-mpa", "mace-ani",
    # ORB 系列
    "orb-v2", "orb-d3-v2", "orb-v3", "orb-v3-mpa", "orb-omat-v3-lora",
    # OMAT24 系列 - 删除这部分
    "omat24-base", "omat24-large", "eqv2-omat", "eqv2-mptrj",
    # GRACE 系列 - 删除这部分
    "grace-2l", "grace-2l-omat", "grace-2m",
    # 其他
    "mattersim-5m",  # 删除此行
    "sevennet-0", "sevennet-mf-ompa", "sevennet-l3i5",  # 删除这些
    "posegnn",  # 删除此行
    "m3gnet", "chgnet"
]

# 修改后（CPU 兼容版本）
CPU_COMPATIBLE_MODELS = [
    # MACE 系列（全部支持）
    "mace-mp", "mace-off23", "mace-omat", "mace-mpa", "mace-ani",
    # ORB 系列（全部支持）
    "orb-v2", "orb-d3-v2", "orb-v3", "orb-v3-mpa", "orb-omat-v3-lora",
    # MatGL 系列
    "m3gnet", "chgnet"
]
```

---

#### 8. `08_microservices_client.py` - 微服务客户端
- ✅ **完全兼容**
- 已支持环境变量配置
- 可以连接任何 MIRA 服务（GPU 或 CPU）

---

## CPU 模式部署

### 测试环境部署

```bash
# 在测试服务器上（如 192.168.100.207）
cd MIRA
./scripts/deploy.sh build-cpu
./scripts/deploy.sh test-cpu
```

**启动的服务：**
- Gateway (端口 8000)
- MACE+ORB Worker (端口 8001) - 支持 MACE 和 ORB 所有模型

### 生产环境部署

```bash
./scripts/deploy.sh up-cpu
```

**启动的服务：**
- Gateway (端口 8000)
- MACE+ORB Worker (端口 8001)
- MatGL Worker (端口 8003) - 支持 M3GNet, CHGNet

---

## 运行 Examples (CPU 模式)

### 方法 1: 本地运行，连接远程 CPU 服务器

```bash
# 设置环境变量指向 CPU 服务器
export MIRA_GATEWAY_URL=http://192.168.100.207:8000

# 运行兼容的示例
python examples/01_basic_usage.py
python examples/03_md_stability.py
python examples/04_bulk_modulus.py
python examples/05_heat_capacity.py
python examples/06_acetylene_adsorption.py
python examples/08_microservices_client.py
```

### 方法 2: 修改后运行需要调整的示例

```bash
# 先手动修改 examples/02_structure_optimization.py
# 将 models 列表改为只包含 CPU 兼容模型

python examples/02_structure_optimization.py

# 或者创建 CPU 版本
cp examples/07_full_benchmark.py examples/07_cpu_benchmark.py
# 编辑 07_cpu_benchmark.py，使用 CPU_COMPATIBLE_MODELS
python examples/07_cpu_benchmark.py
```

### 方法 3: 在服务器上直接运行（推荐）

如果服务器有完整 Python 环境：

```bash
# SSH 到服务器
ssh user@192.168.100.207

# 进入项目目录
cd MIRA

# 确保服务运行中
./scripts/deploy.sh status

# 直接运行（自动连接本地 Gateway）
python examples/01_basic_usage.py
python examples/03_md_stability.py
```

---

## 性能对比

### CPU vs GPU 计算时间对比

以 HKUST-1 (53 atoms) 为例：

| 任务 | GPU 时间 | CPU 时间 | 倍数 |
|------|----------|----------|------|
| 结构优化 (MACE-MP) | ~30s | ~5min | 10x |
| MD 5000 步 | ~2min | ~30min | 15x |
| 体积模量 (9点) | ~1min | ~10min | 10x |
| 声子计算 | ~10min | ~2h | 12x |

### CPU 模式适用场景

✅ **适合：**
- 功能测试和验证
- 小分子计算（<50 原子）
- 单点能量计算
- 开发调试
- 教学演示

❌ **不适合：**
- 大规模基准测试
- 长时间 MD 模拟（>10000 步）
- 大体系计算（>100 原子）
- 生产环境高通量计算

---

## 常见问题

### Q: 为什么某些模型在 CPU 模式下不可用？

A: CPU 模式的 Docker compose 文件只启动了 MACE+ORB 和 MatGL Workers。其他模型（如 SevenNet, GRACE）需要额外的 Workers，目前未包含在 CPU 版本中。

### Q: 能否在 CPU 模式下添加更多模型？

A: 理论上可以，但需要：
1. 创建对应的 CPU Dockerfile
2. 更新 `docker-compose.cpu-prod.yml`
3. 某些模型（如 GRACE）强依赖 GPU 加速，CPU 运行极慢

### Q: CPU 模式下是否所有 examples 都能运行？

A: 只要修改使用 CPU 兼容的模型，所有 examples 都可以运行。关键是选择正确的模型。

### Q: 如何确认服务支持哪些模型？

```bash
# 查询可用模型
curl http://192.168.100.207:8000/api/v1/models | jq
```

---

## 推荐的 CPU 模式工作流

### 1. 开发阶段
```bash
# 本地 CPU 测试
./scripts/deploy.sh test-cpu
python examples/01_basic_usage.py  # 测试连接
python examples/03_md_stability.py  # 功能验证
```

### 2. 远程 CPU 服务器
```bash
# 服务器端
./scripts/deploy.sh up-cpu

# 客户端（Windows/Mac）
export MIRA_GATEWAY_URL=http://your-server:8000
python examples/04_bulk_modulus.py
```

### 3. GPU 生产环境
```bash
# GPU 服务器
./scripts/deploy.sh build
./scripts/deploy.sh up

# 运行完整基准测试
python examples/07_full_benchmark.py
```

---

## 总结

| Example | CPU 兼容 | 需修改 | 说明 |
|---------|----------|--------|------|
| 01_basic_usage.py | ✅ | ❌ | 直接运行 |
| 02_structure_optimization.py | ⚠️ | ✅ | 需修改模型列表 |
| 03_md_stability.py | ✅ | ❌ | 默认 mace-mp，直接运行 |
| 04_bulk_modulus.py | ✅ | ❌ | 默认 mace-mp，直接运行 |
| 05_heat_capacity.py | ✅ | ❌ | 默认 mace-mp，直接运行 |
| 06_acetylene_adsorption.py | ✅ | ❌ | 默认 mace-mp，直接运行 |
| 07_full_benchmark.py | ⚠️ | ✅ | 需大幅修改模型列表 |
| 08_microservices_client.py | ✅ | ❌ | 支持环境变量，直接运行 |

**快速开始（CPU 模式）：**
```bash
# 1. 启动 CPU 服务
./scripts/deploy.sh test-cpu

# 2. 运行兼容示例
for example in 01 03 04 05 06 08; do
    python examples/${example}_*.py
done
```
