# MIRA Examples

本目录包含 MIRA 服务的完整使用示例，覆盖所有功能和模型。

## 目录结构

```
examples/
├── README.md                    # 本文档
├── setup_check.py              # 依赖检查模块
├── config.py                   # 服务器配置模块
├── structures/                  # 示例结构文件
│   ├── HKUST-1.cif             # Cu-BTC MOF
│   ├── MOF-5.cif               # Zn4O(BDC)3
│   ├── ZIF-8.cif               # Zn(mIM)2
│   ├── UiO-66.cif              # Zr-based MOF
│   ├── acetylene.xyz           # C2H2 分子
│   └── HKUST-1_with_C2H2.cif   # 吸附复合物
├── 01_basic_usage.py           # 基础使用
├── 02_structure_optimization.py # 结构优化
├── 03_md_stability.py          # MD 稳定性测试
├── 04_bulk_modulus.py          # 体积模量计算
├── 05_heat_capacity.py         # 热容计算
├── 06_acetylene_adsorption.py  # 乙炔吸附分析
├── 07_full_benchmark.py        # 全模型基准测试
└── 08_microservices_client.py  # 微服务异步客户端
```

## 快速开始

### 1. 安装 ML 力场模型

> ⚠️ **重要**: 不同模型有依赖冲突！建议使用 Docker 微服务或按组合安装。

使用安装脚本（推荐按组合安装）：

```bash
# 检查当前安装状态
python scripts/install_models.py --check

# 推荐组合 A: MACE + ORB (入门首选)
python scripts/install_models.py --combo-a

# 组合 B: FAIRChem + SevenNet
python scripts/install_models.py --combo-b

# 组合 C: MatGL (M3GNet, CHGNet)
python scripts/install_models.py --combo-c

# 组合 D: GRACE
python scripts/install_models.py --combo-d
```

或手动安装：

```bash
pip install mace-torch      # MACE
pip install orb-models      # ORB  
pip install fairchem-core   # FAIRChem/OMAT24
pip install tensorpotential # GRACE
pip install mattersim       # MatterSim
pip install sevenn          # SevenNet
pip install matgl           # MatGL (M3GNet, CHGNet)
```

### 2. 检查依赖

```bash
python examples/setup_check.py
```

### 3. 启动 MIRA 服务

**方式一：本地启动（单环境）**
```bash
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

**方式二：Docker 微服务（推荐，支持所有模型）**
```bash
# GPU 模式
./scripts/deploy.sh build && ./scripts/deploy.sh up

# CPU 模式（无 GPU 环境）
./scripts/deploy.sh build-cpu && ./scripts/deploy.sh up-cpu
```

### 4. 运行示例

```bash
# 本地服务
python examples/01_basic_usage.py

# 连接远程服务器
export MIRA_GATEWAY_URL=http://192.168.100.207:8000
python examples/01_basic_usage.py
```

## 环境变量配置

所有示例都支持通过环境变量配置服务器地址：

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `MIRA_GATEWAY_URL` | `http://localhost:8000` | MIRA Gateway 地址 |

**设置方法：**

```bash
# Linux/macOS
export MIRA_GATEWAY_URL=http://192.168.100.207:8000

# Windows PowerShell
$env:MIRA_GATEWAY_URL = "http://192.168.100.207:8000"
```

## 示例结构

### MOF 材料

| 结构 | 描述 | 特点 |
|-----|------|-----|
| **HKUST-1** | Cu₃(BTC)₂ (Cu-BTC) | 开放金属位点，优异的气体吸附性能 |
| **MOF-5** | Zn₄O(BDC)₃ | 经典 MOF，高比表面积 |
| **ZIF-8** | Zn(mIM)₂ | 沸石咪唑酯框架，高稳定性 |
| **UiO-66** | Zr₆O₄(OH)₄(BDC)₆ | 极高的热/化学/力学稳定性 |

### 客体分子

| 分子 | 描述 | 用途 |
|-----|------|-----|
| **Acetylene (C2H2)** | 乙炔 | 气体分离、吸附研究 |

## 示例说明

### 01. 基础使用 (`01_basic_usage.py`)

演示基本 API 操作：
- 健康检查
- 列出可用模型
- 上传结构文件
- 获取结构信息

```bash
python examples/01_basic_usage.py
```

### 02. 结构优化 (`02_structure_optimization.py`)

使用不同 ML 力场模型进行结构优化：
- 支持 BFGS/FIRE/LBFGS 优化器
- FrechetCellFilter 晶胞优化
- D3 色散校正
- 多模型比较

```bash
python examples/02_structure_optimization.py
```

### 03. MD 稳定性测试 (`03_md_stability.py`)

分子动力学模拟评估结构稳定性：
- NVT 平衡 (Langevin 热浴)
- NPT 生产运行 (NPTBerendsen)
- RMSD 追踪
- 配位数分析
- 温度系列测试

```bash
python examples/03_md_stability.py
```

### 04. 体积模量计算 (`04_bulk_modulus.py`)

E-V 曲线拟合计算力学性质：
- Birch-Murnaghan 状态方程
- 自动 R² 质量评估
- 多模型比较

```bash
python examples/04_bulk_modulus.py
```

### 05. 热容计算 (`05_heat_capacity.py`)

基于声子的热力学性质计算：
- Phonopy 声子计算
- 温度依赖的 Cv
- 虚频检测
- 熵和自由能

```bash
python examples/05_heat_capacity.py
```

### 06. 乙炔吸附分析 (`06_acetylene_adsorption.py`)

MOF-乙炔相互作用能计算：
- 宿主-客体分解
- 能量分解分析
- 多模型比较
- 吸附强度评估

```bash
python examples/06_acetylene_adsorption.py
```

### 07. 全模型基准测试 (`07_full_benchmark.py`)

系统性模型基准测试：
- 20+ 模型全覆盖
- 多结构并行测试
- JSON 结果导出
- 统计汇总报告

```bash
python examples/07_full_benchmark.py
```

## 支持的模型

| 家族 | 模型 | 描述 |
|-----|------|-----|
| **MACE** | mace-mp, mace-off23, mace-omat, mace-mpa, mace-ani | 等变消息传递 |
| **ORB** | orb-v2, orb-d3-v2, orb-v3, orb-v3-mpa, orb-omat-v3-lora | 轨道描述符 |
| **OMAT24** | omat24-base, omat24-large, eqv2-omat, eqv2-mptrj | FAIRChem 系列 |
| **GRACE** | grace-2l, grace-2l-omat, grace-2m | 通用 MLFF |
| **MatterSim** | mattersim-5m | 材料模拟 |
| **SevenNet** | sevennet-0, sevennet-mf-ompa, sevennet-l3i5 | 七体神经网络 |
| **PosEGNN** | posegnn | IBM 位置增强 GNN |
| **MatGL** | m3gnet, chgnet | 图神经网络 |

## 输出示例

### 优化结果
```
=== 模型: mace-mp ===
  最终能量: -1234.5678 eV
  每原子能量: -5.1234 eV/atom
  优化步数: 42
  收敛: True
  RMSD: 0.0234 Å
  体积变化: -1.23%
```

### 稳定性结果
```
模拟时间: 5.00 ps
平均温度: 299.5 K
最终 RMSD: 0.234 Å
体积漂移: 2.1%
稳定性评估: ✓ 结构高度稳定
```

### 体积模量结果
```
体积模量 (B₀): 15.23 GPa
B₀ 导数 (B₀'): 4.21
拟合质量 (R²): 0.999876
```

## 注意事项

1. **GPU 内存**: 大模型（如 MACE-MPA）可能需要 16GB+ 显存
2. **计算时间**: MD 和声子计算可能需要较长时间
3. **结构质量**: 确保输入结构合理，无重叠原子
4. **D3 校正**: 对于 MOF 材料建议启用 D3 色散校正

## 自定义结构

要使用自己的 MOF 结构：

1. 准备 CIF 文件（确保原子坐标正确）
2. 使用 API 上传：
```python
response = requests.post(
    "http://localhost:8000/api/v1/structures/upload",
    data={
        "name": "my-mof",
        "format": "cif",
        "content": open("my-mof.cif").read()
    }
)
structure_id = response.json()['id']
```

## 常见问题

**Q: 服务无法连接？**
- 确保 MIRA 服务已启动
- 检查端口是否正确 (默认 8000)
- 如连接远程服务器，检查环境变量 `MIRA_GATEWAY_URL` 是否设置正确

**Q: 模型加载失败？**
- 检查 GPU 显存是否充足
- 确保相关模型包已安装
- 使用 `python scripts/install_models.py --check` 检查安装状态

**Q: 计算超时？**
- 增加 timeout 参数
- 减少 MD 步数或优化步数

**Q: 依赖冲突？**
- 使用 Docker 微服务架构: `./scripts/deploy.sh up`
- 或使用组合安装: `--combo-a`, `--combo-b` 等

## 联系

**作者**: 李世博 (Shibo Li)  
**邮箱**: shadow.li981@gmail.com
