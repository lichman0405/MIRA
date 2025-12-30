# MIRA 微服务架构部署指南

## 架构概述

MIRA 采用微服务架构，解决 ML 力场模型之间的依赖冲突问题。

```
┌─────────────────────────────────────────────────────────────┐
│                      用户请求                                │
│                   http://localhost:8000                      │
└─────────────────────────┬───────────────────────────────────┘
                          ▼
┌─────────────────────────────────────────────────────────────┐
│                 MIRA Gateway (FastAPI)                       │
│                    主服务 - 端口 8000                         │
│  • 接收用户请求                                               │
│  • 自动路由到对应模型服务                                      │
│  • 聚合结果返回                                               │
│  • 支持多模型并行计算                                          │
└──────┬──────────┬──────────┬──────────┬──────────┬──────────┘
       │          │          │          │          │
       ▼          ▼          ▼          ▼          ▼
┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
│ MACE+ORB │ │FAIRChem+ │ │  MatGL   │ │  GRACE   │ │MatterSim │
│   :8001  │ │SevenNet  │ │  :8003   │ │  :8004   │ │  :8005   │
│          │ │  :8002   │ │          │ │          │ │          │
│ GPU 0    │ │  GPU 1   │ │  GPU 2   │ │  GPU 3   │ │  GPU 4   │
│ e3nn=0.4 │ │ e3nn>=0.5│ │   DGL    │ │TensorFlow│ │ PyTorch  │
└──────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘
   容器 A        容器 B       容器 C       容器 D       容器 E
```

## 快速开始

### 1. GPU 测试环境部署 (单 GPU)

```bash
cd /path/to/MIRA

# 构建镜像
./scripts/deploy.sh build

# 启动测试环境 (只有 Gateway + MACE-ORB)
./scripts/deploy.sh test

# 检查状态
./scripts/deploy.sh status

# 访问 API 文档
# http://localhost:8000/docs
```

### 2. GPU 生产环境部署 (多 GPU)

```bash
# 构建所有镜像
./scripts/deploy.sh build

# 启动所有服务
./scripts/deploy.sh up

# 查看日志
./scripts/deploy.sh logs

# 停止服务
./scripts/deploy.sh down
```

### 3. CPU 模式部署 (无 GPU 环境)

适用于没有 NVIDIA GPU 的测试服务器或开发环境。

```bash
# 构建 CPU 镜像
./scripts/deploy.sh build-cpu

# 启动 CPU 测试环境 (Gateway + MACE-ORB)
./scripts/deploy.sh test-cpu

# 启动 CPU 生产环境 (Gateway + MACE-ORB + MatGL)
./scripts/deploy.sh up-cpu

# 检查状态
./scripts/deploy.sh status
```

**CPU 模式注意事项：**
- ⚠️ 计算速度比 GPU 慢 10-100 倍
- ⚠️ 适合功能验证和小批量计算
- ⚠️ 目前支持 MACE、ORB、MatGL 模型
- ✅ 无需 NVIDIA 驱动和 CUDA

### 4. 手动 Docker Compose

```bash
cd docker

# GPU 测试环境
docker compose -f docker-compose.test.yml up -d

# GPU 生产环境 (需要 5+ GPU)
docker compose -f docker-compose.microservices.yml up -d

# CPU 测试环境
docker compose -f docker-compose.cpu.yml up -d

# CPU 生产环境
docker compose -f docker-compose.cpu-prod.yml up -d
```

## GPU 分配

默认 GPU 分配策略（8 卡服务器）:

| GPU | 服务 | 模型 |
|-----|------|------|
| 0 | mace-orb | MACE, ORB |
| 1 | fairchem-sevennet | OMAT24, SevenNet |
| 2 | matgl | M3GNet, CHGNet |
| 3 | grace | GRACE |
| 4 | mattersim | MatterSim |
| 5-7 | 备用/扩展 | - |

### 自定义 GPU 分配

编辑 `docker/docker-compose.microservices.yml`:

```yaml
mace-orb:
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            device_ids: ['0']  # 修改为需要的 GPU ID
            capabilities: [gpu]
```

## API 使用

### 健康检查

```bash
curl http://localhost:8000/health
```

### 列出可用模型

```bash
curl http://localhost:8000/api/v1/models
```

### 单点能量计算

```bash
curl -X POST http://localhost:8000/api/v1/single_point \
  -H "Content-Type: application/json" \
  -d '{
    "model_name": "mace-mp",
    "atoms": {
      "symbols": ["Cu", "Cu", "Cu", "Cu"],
      "positions": [[0,0,0], [1.8,1.8,0], [1.8,0,1.8], [0,1.8,1.8]],
      "cell": [[3.6,0,0], [0,3.6,0], [0,0,3.6]],
      "pbc": [true, true, true]
    }
  }'
```

### 多模型比较

```bash
curl -X POST http://localhost:8000/api/v1/multi_model \
  -H "Content-Type: application/json" \
  -d '{
    "model_names": ["mace-mp", "m3gnet", "chgnet"],
    "task": "single_point",
    "atoms": {
      "symbols": ["Cu", "Cu", "Cu", "Cu"],
      "positions": [[0,0,0], [1.8,1.8,0], [1.8,0,1.8], [0,1.8,1.8]],
      "cell": [[3.6,0,0], [0,3.6,0], [0,0,3.6]],
      "pbc": [true, true, true]
    }
  }'
```

## 服务管理

### 查看日志

```bash
# 所有服务
./scripts/deploy.sh logs

# 特定服务
./scripts/deploy.sh logs gateway
./scripts/deploy.sh logs mace-orb
```

### 重启单个服务

```bash
cd docker
docker compose -f docker-compose.microservices.yml restart mace-orb
```

### 扩展服务副本

```bash
cd docker
docker compose -f docker-compose.microservices.yml up -d --scale mace-orb=2
```

## 故障排查

### 1. 服务无法启动

```bash
# 查看详细日志
docker logs mira-mace-orb

# 检查 GPU 可用性
nvidia-smi
docker run --rm --gpus all nvidia/cuda:12.1-base nvidia-smi
```

### 2. Gateway 无法连接 Worker

```bash
# 检查网络
docker network inspect docker_mira-network

# 检查服务是否在同一网络
docker inspect mira-gateway | grep NetworkMode
docker inspect mira-mace-orb | grep NetworkMode
```

### 3. GPU 内存不足

```bash
# 查看 GPU 使用情况
nvidia-smi

# 清理未使用的容器
docker system prune -f
```

## 文件结构

```
MIRA/
├── docker/
│   ├── Dockerfile.base              # 基础镜像
│   ├── Dockerfile.gateway           # Gateway 镜像
│   ├── Dockerfile.mace-orb          # MACE+ORB GPU 镜像
│   ├── Dockerfile.mace-orb-cpu      # MACE+ORB CPU 镜像
│   ├── Dockerfile.fairchem-sevennet # FAIRChem+SevenNet 镜像
│   ├── Dockerfile.matgl             # MatGL GPU 镜像
│   ├── Dockerfile.matgl-cpu         # MatGL CPU 镜像
│   ├── Dockerfile.grace             # GRACE 镜像
│   ├── Dockerfile.mattersim         # MatterSim 镜像
│   ├── docker-compose.microservices.yml  # GPU 生产环境
│   ├── docker-compose.test.yml      # GPU 测试环境
│   ├── docker-compose.cpu.yml       # CPU 测试环境
│   └── docker-compose.cpu-prod.yml  # CPU 生产环境
├── services/
│   ├── shared/                      # 共享代码
│   │   ├── schemas.py               # 统一数据模型
│   │   └── worker_base.py           # Worker 基类
│   ├── gateway/                     # Gateway 服务
│   │   └── main.py
│   ├── mace_orb/                    # MACE+ORB Worker
│   │   └── main.py
│   ├── fairchem_sevennet/           # FAIRChem+SevenNet Worker
│   │   └── main.py
│   ├── matgl/                       # MatGL Worker
│   │   └── main.py
│   ├── grace/                       # GRACE Worker
│   │   └── main.py
│   └── mattersim/                   # MatterSim Worker
│       └── main.py
└── scripts/
    └── deploy.sh                    # 部署脚本
```

## 部署脚本命令参考

| 命令 | 说明 |
|------|------|
| `./scripts/deploy.sh build` | 构建所有 GPU 镜像 |
| `./scripts/deploy.sh build-cpu` | 构建 CPU 镜像 |
| `./scripts/deploy.sh test` | 启动 GPU 测试环境 |
| `./scripts/deploy.sh test-cpu` | 启动 CPU 测试环境 |
| `./scripts/deploy.sh up` | 启动 GPU 生产环境 |
| `./scripts/deploy.sh up-cpu` | 启动 CPU 生产环境 |
| `./scripts/deploy.sh down` | 停止所有服务 |
| `./scripts/deploy.sh logs` | 查看日志 |
| `./scripts/deploy.sh logs <服务名>` | 查看特定服务日志 |
| `./scripts/deploy.sh status` | 查看服务状态 |
| `./scripts/deploy.sh clean` | 清理 Docker 资源 |

## 注意事项

1. **首次构建耗时较长**: 每个模型容器需要安装各自的依赖，首次构建可能需要 30-60 分钟
2. **模型下载**: 首次运行时，模型会自动下载到 `mira-models` 卷
3. **GPU 内存**: 每个模型服务约需 4-8GB GPU 内存
4. **网络**: 所有服务通过 Docker 内部网络通信，Gateway 对外暴露 8000 端口
5. **CPU 模式**: 计算速度较慢，适合功能测试和小批量计算

## 客户端使用

### 远程连接

客户端机器无需安装 ML 模型包，只需 `requests` 和 `ase` 库：

```bash
pip install requests ase
```

### 环境变量配置

```bash
# Linux/macOS
export MIRA_GATEWAY_URL=http://192.168.100.207:8000

# Windows PowerShell
$env:MIRA_GATEWAY_URL = "http://192.168.100.207:8000"

# 或在 Python 中设置
import os
os.environ["MIRA_GATEWAY_URL"] = "http://192.168.100.207:8000"
```

### 运行示例

```bash
cd MIRA/examples

# 基础使用
python 01_basic_usage.py

# 结构优化
python 02_structure_optimization.py

# 多模型基准测试
python 07_multi_model_benchmark.py
```

### 部署配置对比

| 配置 | 命令 | 服务 | 适用场景 |
|------|------|------|----------|
| GPU 测试 | `deploy.sh test` | Gateway + MACE-ORB | 开发调试 |
| GPU 生产 | `deploy.sh up` | 全部 5 个 Worker | 正式计算 |
| CPU 测试 | `deploy.sh test-cpu` | Gateway + MACE-ORB (CPU) | 无 GPU 功能测试 |
| CPU 生产 | `deploy.sh up-cpu` | Gateway + MACE-ORB + MatGL (CPU) | 无 GPU 小批量计算 |
