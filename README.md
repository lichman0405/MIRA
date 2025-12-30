# MIRA - MiQroEra Interatomic-potential Reliability Arena

一个用于在金属有机框架 (MOFs) 上对机器学习原子间势进行全面基准测试的 FastAPI 服务。

## 🏗️ 系统架构

MIRA 采用 **微服务架构**，每个模型家族运行在独立的 Docker 容器中，完美解决依赖冲突问题。

```
┌─────────────────────────────────────────────────────────────┐
│                      用户请求                                │
│                   http://localhost:8000                      │
└─────────────────────────┬───────────────────────────────────┘
                          ▼
┌─────────────────────────────────────────────────────────────┐
│                 MIRA Gateway (FastAPI)                       │
│                    主服务 - 端口 8000                         │
│  • 统一入口，自动路由                                         │
│  • 支持多模型并行计算                                         │
└──────┬──────────┬──────────┬──────────┬──────────┬──────────┘
       │          │          │          │          │
       ▼          ▼          ▼          ▼          ▼
┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
│ MACE+ORB │ │FAIRChem+ │ │  MatGL   │ │  GRACE   │ │MatterSim │
│   :8001  │ │SevenNet  │ │  :8003   │ │  :8004   │ │  :8005   │
│  GPU 0   │ │  :8002   │ │  GPU 2   │ │  GPU 3   │ │  GPU 4   │
│ e3nn=0.4 │ │  GPU 1   │ │   DGL    │ │TensorFlow│ │ PyTorch  │
└──────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘
```

## ✨ 功能特性

### 支持的模型 (20+ 变体)

| 模型家族 | 模型变体 | 描述 |
|----------|----------|------|
| **MACE** | MP, OFF23, OMAT, MPA, ANI | 等变消息传递网络 |
| **ORB** | v2, v3, OMAT-v3-LoRA | 轨道基描述符 |
| **OMAT24** | OMat, eqV2 变体 | Meta FAIRChem 模型 |
| **GRACE** | 2L, 2M | 通用 ML 力场 |
| **MatterSim** | 5M | 材料模拟 |
| **SevenNet** | 0, MF-ompa, l3i5 | 七体神经网络 |
| **PosEGNN** | IBM model | 位置增强 GNN |
| **MatGL** | M3GNet, CHGNet | 图神经网络 |

### 计算任务

1. **结构优化**
   - 优化器: BFGS, FIRE, LBFGS
   - 晶胞过滤器: FrechetCellFilter 用于全松弛
   - 支持 D3 色散校正

2. **MD 稳定性测试**
   - NVT 平衡 (Langevin 恒温器)
   - NPT 生产 (NPTBerendsen)
   - 配位数分析
   - RMSD 追踪

3. **体积模量计算**
   - E-V 曲线采样
   - Birch-Murnaghan 状态方程拟合
   - 自动 R² 质量评估

4. **热容计算**
   - 使用 Phonopy 进行声子计算
   - 温度依赖的 Cv
   - 虚频模式检测

5. **QMOF 能量评估**
   - 单点能量计算
   - 与 DFT 参考值比较

6. **相互作用能分析**
   - 主客体分解
   - 组分能量分解

## 📦 安装

### 环境要求

- Python >= 3.10
- CUDA >= 12.0 (用于 GPU 加速)
- ASE >= 3.27.0 (包含 NPT 支持)

### 🐳 Docker 微服务部署 (推荐)

```bash
# 克隆仓库
git clone https://github.com/lichman0405/MIRA.git
cd MIRA

# 构建 Docker 镜像
./scripts/deploy.sh build        # GPU 版本
./scripts/deploy.sh build-cpu    # CPU 版本 (无 GPU 环境)

# 启动测试环境 (单 GPU: Gateway + MACE-ORB)
./scripts/deploy.sh test         # GPU 模式
./scripts/deploy.sh test-cpu     # CPU 模式

# 启动生产环境 (多 GPU: 所有服务)
./scripts/deploy.sh up

# 查看状态
./scripts/deploy.sh status

# 访问 API 文档
# http://localhost:8000/docs
```

详细部署指南: [docs/DEPLOYMENT.md](docs/DEPLOYMENT.md)

### 快速开始 (传统方式)

```bash
# 克隆仓库
git clone https://github.com/lichman0405/MIRA.git
cd MIRA

# 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Linux/macOS
# 或: venv\Scripts\activate  # Windows

# 安装基础依赖
pip install -r requirements.txt

# 安装 ML 力场模型
python scripts/install_models.py --check      # 检查状态
python scripts/install_models.py --combo-a    # MACE + ORB (推荐)

# 运行服务
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
```

### ML 力场安装

> ⚠️ **重要提示**: 不同 ML 力场模型有不兼容的依赖版本！建议使用 Docker 微服务或多 conda 环境策略。

**兼容的模型组合：**

| 组合 | 模型 | 依赖 | 适用场景 |
|------|------|------|----------|
| 🅰 A | MACE + ORB | PyTorch + e3nn==0.4.4 | 推荐入门、MOF 基准测试 |
| 🅱 B | FAIRChem + SevenNet | PyTorch + e3nn>=0.5 | 大规模材料预测 |
| 🅲 C | MatGL | PyTorch + DGL | 电池材料、晶体结构 |
| 🅳 D | GRACE | TensorFlow | 高精度力场 |

**推荐安装方式（多环境）：**

```bash
# 环境 1: MACE + ORB (推荐入门)
conda create -n mira-mace python=3.10
conda activate mira-mace
pip install -e .
python scripts/install_models.py --combo-a

# 环境 2: FAIRChem + SevenNet
conda create -n mira-fairchem python=3.10
conda activate mira-fairchem
pip install -e .
python scripts/install_models.py --combo-b
```

**使用安装脚本：**

```bash
# 检查已安装的模型
python scripts/install_models.py --check

# 安装推荐组合 (MACE + ORB)
python scripts/install_models.py --combo-a

# 安装其他组合
python scripts/install_models.py --combo-b  # FAIRChem + SevenNet
python scripts/install_models.py --combo-c  # MatGL
python scripts/install_models.py --combo-d  # GRACE

# 安装单个模型
python scripts/install_models.py --mace
python scripts/install_models.py --mace --orb
```

**手动安装单个模型：**

```bash
pip install mace-torch      # MACE
pip install orb-models      # ORB  
pip install fairchem-core   # FAIRChem/OMAT24
pip install tensorpotential # GRACE
pip install mattersim       # MatterSim
pip install sevenn          # SevenNet
pip install matgl           # MatGL (M3GNet, CHGNet)
```

## 📖 API 使用

### 基础 URL
```
http://localhost:8000/api/v1
```

### 交互式文档
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

### 环境变量配置

连接远程服务器时，设置环境变量：

```bash
# Linux/macOS
export MIRA_GATEWAY_URL=http://192.168.100.207:8000

# Windows PowerShell
$env:MIRA_GATEWAY_URL = "http://192.168.100.207:8000"
```

### 示例：结构优化

```python
import os
import requests

# 获取服务地址
BASE_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000") + "/api/v1"

# 上传结构
with open("structure.cif", "r") as f:
    content = f.read()

response = requests.post(
    f"{BASE_URL}/structures/upload",
    data={
        "name": "MOF-5",
        "format": "cif",
        "content": content
    }
)
structure_id = response.json()["id"]

# 提交优化任务
response = requests.post(
    f"{BASE_URL}/tasks/optimization",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "fmax": 0.01,
        "max_steps": 500,
        "optimizer": "BFGS",
        "use_filter": True,
        "enable_d3": True
    }
)
task_id = response.json()["task_id"]

# 检查进度
response = requests.get(f"{BASE_URL}/tasks/{task_id}")
print(response.json())

# 获取结果
response = requests.get(f"{BASE_URL}/results/{task_id}")
result = response.json()
print(f"最终能量: {result['final_energy']} eV")
```

### 示例：MD 稳定性测试

```python
response = requests.post(
    f"{BASE_URL}/tasks/stability",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "nvt_steps": 5000,
        "nvt_temperature": 300,
        "npt_steps": 10000,
        "npt_temperature": 300,
        "npt_pressure": 1.01325,
        "timestep": 1.0,
        "pre_optimize": True,
        "enable_d3": True
    }
)
```

## ⚙️ 配置

环境变量（或 `.env` 文件）：

| 变量名 | 默认值 | 说明 |
|--------|--------|------|
| `DEBUG` | false | 启用调试模式 |
| `DEFAULT_DEVICE` | cuda | 默认计算设备 |
| `STRUCTURES_DIR` | ./data/structures | 结构存储目录 |
| `RESULTS_DIR` | ./data/results | 结果存储目录 |
| `MAX_WORKERS` | 4 | 并行任务工作线程数 |
| `MAX_MD_STEPS` | 100000 | MD 步数限制 |
| `MAX_OPT_STEPS` | 2000 | 优化步数限制 |
| `CORS_ORIGINS` | ["*"] | 允许的 CORS 来源 |

## 📁 项目结构

```
MIRA/
├── app/                      # 主应用代码
│   ├── main.py              # FastAPI 应用入口
│   ├── config.py            # 配置管理
│   ├── api/v1/              # API 路由
│   ├── schemas/             # 数据模式
│   ├── models/              # 模型适配器
│   ├── services/            # 业务服务
│   └── core/                # 核心工具
├── services/                 # 微服务工作节点
│   ├── shared/              # 共享代码
│   ├── gateway/             # API 网关
│   ├── mace_orb/            # MACE+ORB Worker
│   ├── fairchem_sevennet/   # FAIRChem+SevenNet Worker
│   ├── matgl/               # MatGL Worker
│   ├── grace/               # GRACE Worker
│   └── mattersim/           # MatterSim Worker
├── docker/                   # Docker 配置
│   ├── Dockerfile.*         # 各服务镜像
│   ├── docker-compose.microservices.yml  # GPU 生产
│   ├── docker-compose.test.yml           # GPU 测试
│   ├── docker-compose.cpu.yml            # CPU 测试
│   └── docker-compose.cpu-prod.yml       # CPU 生产
├── scripts/
│   ├── install_models.py    # ML 力场安装脚本
│   └── deploy.sh            # 部署脚本
├── examples/
│   ├── 01-07_*.py           # 功能示例
│   ├── 08_microservices_client.py  # 异步客户端
│   ├── config.py            # 服务器配置
│   └── structures/          # 示例 MOF 结构
├── docs/
│   └── DEPLOYMENT.md        # 部署文档
├── .env.example             # 环境变量模板
├── requirements.txt         # Python 依赖
└── README.md
```

## 📚 示例

查看 [examples/](examples/) 目录获取完整使用示例：

```bash
# 首先检查依赖
python examples/setup_check.py

# 运行示例
python examples/01_basic_usage.py
python examples/02_structure_optimization.py
python examples/07_full_benchmark.py

# 连接远程服务器运行
export MIRA_GATEWAY_URL=http://192.168.100.207:8000
python examples/01_basic_usage.py
```

## 📝 注意事项

- **GPU 内存**: 大型模型 (MACE-MPA, SevenNet-l3i5) 可能需要 16GB+ 显存
- **ASE 版本**: 需要 ASE >= 3.27.0 以支持 NPT 动力学
- **D3 校正**: 某些任务可从 DFT-D3 色散校正中受益
- **声子计算**: 热容计算需要足够的超胞尺寸

## 👤 作者

**李世博** (Shibo Li)  
📧 shadow.li981@gmail.com

## 📄 许可证

MIT License

## 📖 引用

如果您在研究中使用了 MIRA，请引用：

```bibtex
@software{mira2025,
  author={Li, Shibo},
  title={MIRA: MiQroEra Interatomic-potential Reliability Arena},
  year={2025},
  url={https://github.com/lichman0405/MIRA}
}
```
