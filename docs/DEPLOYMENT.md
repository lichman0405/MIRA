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

### 1. 测试环境部署 (单 GPU)

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

### 2. 生产环境部署 (多 GPU)

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

### 3. 手动 Docker Compose

```bash
cd docker

# 测试环境
docker compose -f docker-compose.test.yml up -d

# 生产环境 (需要 5+ GPU)
docker compose -f docker-compose.microservices.yml up -d
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
│   ├── Dockerfile.gateway           # Gateway 镜像
│   ├── Dockerfile.mace-orb          # MACE+ORB 镜像
│   ├── Dockerfile.fairchem-sevennet # FAIRChem+SevenNet 镜像
│   ├── Dockerfile.matgl             # MatGL 镜像
│   ├── Dockerfile.grace             # GRACE 镜像
│   ├── Dockerfile.mattersim         # MatterSim 镜像
│   ├── docker-compose.microservices.yml  # 生产环境
│   └── docker-compose.test.yml      # 测试环境
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

## 注意事项

1. **首次构建耗时较长**: 每个模型容器需要安装各自的依赖，首次构建可能需要 30-60 分钟
2. **模型下载**: 首次运行时，模型会自动下载到 `mira-models` 卷
3. **GPU 内存**: 每个模型服务约需 4-8GB GPU 内存
4. **网络**: 所有服务通过 Docker 内部网络通信，Gateway 对外暴露 8000 端口
