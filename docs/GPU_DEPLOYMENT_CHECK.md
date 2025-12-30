# MIRA 多 GPU 生产环境配置检查报告

## 执行时间
2025-12-30

## 检查范围
- Docker Compose 配置
- GPU 分配策略
- Dockerfile 配置
- 部署脚本
- README 文档一致性

---

## ✅ 检查结果总结

### 1. GPU 分配配置 ✅ 正确

**docker-compose.microservices.yml 中的 GPU 分配：**

| 服务 | 容器名 | 端口 | GPU ID | 环境变量 | 状态 |
|------|--------|------|--------|----------|------|
| gateway | mira-gateway | 8000 | - | - | ✅ 无 GPU 需求 |
| mace-orb | mira-mace-orb | 8001 | 0 | `NVIDIA_VISIBLE_DEVICES=0` | ✅ 正确 |
| fairchem-sevennet | mira-fairchem-sevennet | 8002 | 1 | `NVIDIA_VISIBLE_DEVICES=1` | ✅ 正确 |
| matgl | mira-matgl | 8003 | 2 | `NVIDIA_VISIBLE_DEVICES=2` | ✅ 正确 |
| grace | mira-grace | 8004 | 3 | `NVIDIA_VISIBLE_DEVICES=3` | ✅ 正确 |
| mattersim | mira-mattersim | 8005 | 4 | `NVIDIA_VISIBLE_DEVICES=4` | ✅ 正确 |

**配置要点：**
- ✅ 每个 Worker 独立 GPU（避免资源竞争）
- ✅ `device_ids` 和 `NVIDIA_VISIBLE_DEVICES` 一致
- ✅ 容器内部都映射为 `CUDA_VISIBLE_DEVICES=0`（标准做法）
- ✅ 使用 5 个 GPU（GPU 0-4），适合 8 卡服务器

### 2. 与 README 的一致性检查 ⚠️ 需要更新

**README.md 中的架构图：**
```
┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
│ MACE+ORB │ │FAIRChem+ │ │  MatGL   │ │  GRACE   │ │MatterSim │
│   :8001  │ │SevenNet  │ │  :8003   │ │  :8004   │ │  :8005   │
│  GPU 0   │ │  :8002   │ │  GPU 2   │ │  GPU 3   │ │  GPU 4   │
│ e3nn=0.4 │ │  GPU 1   │ │   DGL    │ │TensorFlow│ │ PyTorch  │
└──────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘
```

**实际配置：**
- ✅ MACE+ORB: GPU 0, 端口 8001
- ✅ FAIRChem+SevenNet: GPU 1, 端口 8002
- ✅ MatGL: GPU 2, 端口 8003
- ✅ GRACE: GPU 3, 端口 8004
- ✅ MatterSim: GPU 4, 端口 8005

**结论：README 架构图与实际配置 100% 一致！** ✅

### 3. Dockerfile 配置检查 ✅ 正确

所有 GPU Dockerfile 都：
- ✅ 使用 `nvidia/cuda:12.1-cudnn8-runtime-ubuntu22.04` 基础镜像
- ✅ 安装 PyTorch with CUDA 12.1
- ✅ ASE >= 3.27.0（已修复导入问题）
- ✅ 包含健康检查
- ✅ 正确的端口映射
- ✅ 复制共享模块和服务代码

### 4. 部署脚本检查 ✅ 正确

**scripts/deploy.sh 支持的命令：**

| 命令 | 功能 | 文件 | 状态 |
|------|------|------|------|
| `build` | 构建所有 GPU 镜像 | - | ✅ 正确 |
| `build-cpu` | 构建 CPU 镜像 | - | ✅ 正确 |
| `test` | GPU 测试环境 | docker-compose.test.yml | ✅ 正确 |
| `test-cpu` | CPU 测试环境 | docker-compose.cpu.yml | ✅ 正确 |
| `up` | GPU 生产环境 | docker-compose.microservices.yml | ✅ 正确 |
| `up-cpu` | CPU 生产环境 | docker-compose.cpu-prod.yml | ✅ 正确 |
| `down` | 停止所有服务 | 所有 compose 文件 | ✅ 正确 |
| `status` | 查看状态 | - | ✅ 正确 |
| `logs` | 查看日志 | - | ✅ 正确 |
| `clean` | 清理资源 | - | ✅ 正确 |

### 5. 网络和存储配置 ✅ 正确

**网络：**
```yaml
networks:
  mira-network:
    driver: bridge
```
- ✅ 所有服务在同一网络
- ✅ 服务间通信使用容器名（如 `mira-mace-orb:8001`）

**存储卷：**
```yaml
volumes:
  mira-models:
    driver: local
```
- ✅ 共享模型存储卷
- ✅ 避免每个容器重复下载模型

### 6. 健康检查配置 ✅ 完善

所有服务都配置了健康检查：
- ✅ Gateway: 30s 间隔，10s 启动期
- ✅ Workers: 30s 间隔，120s 启动期（考虑模型加载时间）
- ✅ 3 次重试
- ✅ curl 检查 `/health` 端点

### 7. 日志配置 ✅ 合理

```yaml
logging:
  driver: "json-file"
  options:
    max-size: "100m"
    max-file: "3"
```
- ✅ 限制日志文件大小（100MB）
- ✅ 保留 3 个历史文件
- ✅ 防止磁盘空间耗尽

---

## ⚠️ 发现的潜在问题

### 1. GPU 内存分配策略 (建议优化)

当前每个 Worker 独占一个 GPU，这是**安全但可能不够高效**的策略。

**建议优化方案：**

#### 方案 A: 基于模型大小的 GPU 共享（推荐）

某些轻量级模型可以共享 GPU：

```yaml
# 轻量级模型可以共享 GPU
mace-orb:
  device_ids: ['0']  # MACE-MP, ORB-v2 较小

matgl:
  device_ids: ['0']  # M3GNet, CHGNet 也较小，可与 MACE 共享

fairchem-sevennet:
  device_ids: ['1', '2']  # OMAT24, SevenNet 较大，独占或双卡

grace:
  device_ids: ['3']  # GRACE 独占

mattersim:
  device_ids: ['4']  # MatterSim 独占
```

#### 方案 B: 保持当前配置（最稳定）

如果不确定内存需求，**当前配置是最安全的选择**。

### 2. docker-compose 文件位置

**当前位置：** `docker/docker-compose.*.yml`  
**部署脚本引用：** `$DOCKER_DIR/docker-compose.*.yml`

✅ **一致，无问题**

但是，常见的最佳实践是将 compose 文件放在项目根目录：

```
MIRA/
├── docker-compose.yml          # 主配置
├── docker-compose.prod.yml     # 生产环境覆盖
├── docker-compose.test.yml     # 测试环境
└── docker/
    └── Dockerfile.*            # 只放 Dockerfile
```

**建议：** 当前配置也可以，如果要修改需要同步更新部署脚本。

### 3. 缺少环境变量文件

建议创建 `.env.production` 用于生产环境配置：

```bash
# .env.production
COMPOSE_PROJECT_NAME=mira
MIRA_ENV=production
CUDA_VERSION=12.1
```

然后在 compose 文件中引用：

```yaml
env_file:
  - ../.env.production
```

---

## 🎯 针对 8 GPU 服务器的推荐配置

假设您有 8 张 GPU（如 8x A100），以下是优化建议：

### 配置 1: 保守策略（当前配置）✅

```
GPU 0: MACE-ORB
GPU 1: FAIRChem-SevenNet
GPU 2: MatGL
GPU 3: GRACE
GPU 4: MatterSim
GPU 5-7: 保留/扩展
```

**优点：**
- ✅ 最稳定，无资源竞争
- ✅ 故障隔离性好
- ✅ 易于调试

**缺点：**
- ⚠️ GPU 利用率可能不高（小模型浪费）
- ⚠️ 3 张 GPU 闲置

### 配置 2: 高效策略（推荐）

```yaml
# 轻量级服务共享 GPU
GPU 0: MACE-ORB (单卡足够)
GPU 1: MatGL (单卡足够)

# 中等服务独占 GPU
GPU 2: FAIRChem Worker 1
GPU 3: SevenNet Worker 1
GPU 4: GRACE
GPU 5: MatterSim

# 负载均衡副本（高并发场景）
GPU 6: FAIRChem Worker 2 (副本)
GPU 7: SevenNet Worker 2 (副本)
```

**实现方式：**

```yaml
# 添加副本服务
fairchem-sevennet-replica:
  <<: *fairchem-sevennet  # 继承配置
  container_name: mira-fairchem-sevennet-2
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            device_ids: ['6']
            capabilities: [gpu]
  environment:
    - NVIDIA_VISIBLE_DEVICES=6
```

### 配置 3: 动态扩展策略（最灵活）

使用 Docker Compose 的 `scale` 功能：

```bash
# 根据负载动态扩展
docker-compose up -d --scale mace-orb=2
docker-compose up -d --scale fairchem-sevennet=3
```

但需要配置负载均衡器（如 Nginx）。

---

## 🔍 测试建议

### 1. GPU 可用性测试

```bash
# 在每个 Worker 容器中测试
docker exec mira-mace-orb python -c "import torch; print(f'GPU 可用: {torch.cuda.is_available()}, 设备数: {torch.cuda.device_count()}')"
docker exec mira-fairchem-sevennet python -c "import torch; print(f'GPU: {torch.cuda.get_device_name(0)}')"
```

### 2. GPU 内存使用监控

```bash
# 在宿主机监控
watch -n 1 nvidia-smi

# 或使用 GPUtil
python -c "import GPUtil; GPUtil.showUtilization()"
```

### 3. 并发负载测试

```bash
# 并发提交 10 个优化任务
for i in {1..10}; do
  curl -X POST http://localhost:8000/api/v1/tasks/optimization \
    -H "Content-Type: application/json" \
    -d @test_task.json &
done
```

### 4. Worker 故障恢复测试

```bash
# 模拟 Worker 崩溃
docker stop mira-mace-orb

# 检查 Gateway 是否正确处理
curl http://localhost:8000/health

# 重启 Worker
docker start mira-mace-orb
```

---

## 📋 部署检查清单

启动生产环境前检查：

- [ ] ✅ Docker 和 Docker Compose 已安装
- [ ] ✅ NVIDIA Docker Runtime 已配置
- [ ] ✅ GPU 驱动和 CUDA 12.1+ 已安装
- [ ] ✅ 运行 `nvidia-smi` 确认所有 GPU 可见
- [ ] ✅ 拉取最新代码（包含 ASE 修复）
- [ ] ✅ 构建镜像：`./scripts/deploy.sh build`
- [ ] ✅ 启动服务：`./scripts/deploy.sh up`
- [ ] ✅ 检查健康状态：`./scripts/deploy.sh status`
- [ ] ✅ 查看日志确认无错误
- [ ] ✅ 测试 API：`curl http://localhost:8000/api/v1/models`
- [ ] ✅ 运行 examples 验证功能

---

## ✅ 最终结论

### 当前配置评估：**A 级（优秀）**

**优点：**
1. ✅ GPU 分配策略清晰合理
2. ✅ 架构图与实际配置完全一致
3. ✅ Dockerfile 配置规范
4. ✅ 部署脚本功能完善
5. ✅ 健康检查和日志配置完善
6. ✅ 网络和存储设计合理
7. ✅ ASE 兼容性问题已修复

**建议改进（可选）：**
1. ⚡ 考虑 GPU 共享策略（提高利用率）
2. 📊 添加监控（Prometheus + Grafana）
3. 🔄 配置自动重启策略（已有 `restart: unless-stopped`）
4. 📝 添加 `.env.production` 环境配置文件

### 对于 8 GPU 生产服务器

**当前配置完全可用！** 建议：

1. **初期部署**：使用当前配置（最稳定）
2. **监控阶段**：观察 GPU 利用率
3. **优化阶段**：根据实际负载调整 GPU 分配

**快速启动命令：**

```bash
cd MIRA
git pull                          # 拉取最新修复
./scripts/deploy.sh build         # 构建镜像
./scripts/deploy.sh up            # 启动所有服务
./scripts/deploy.sh status        # 检查状态
```

**预期结果：**
- ✅ 5 个 Worker 容器各占用一张 GPU (GPU 0-4)
- ✅ Gateway 不使用 GPU
- ✅ GPU 5-7 空闲（可用于未来扩展）
- ✅ 所有服务健康，API 可访问
- ✅ 无 ASE 导入错误
- ✅ 支持所有 20+ 模型变体

---

## 补充：完整的生产部署流程

```bash
# 1. 准备环境
nvidia-smi  # 确认 8 张 GPU 都可见

# 2. 克隆或更新代码
cd /path/to/MIRA
git pull

# 3. 构建镜像（约 30-60 分钟）
./scripts/deploy.sh build

# 4. 启动服务
./scripts/deploy.sh up

# 5. 等待所有服务就绪（约 2-5 分钟）
watch -n 5 ./scripts/deploy.sh status

# 6. 验证 API
curl http://localhost:8000/health | jq
curl http://localhost:8000/api/v1/models | jq

# 7. 运行测试
cd examples
export MIRA_GATEWAY_URL=http://localhost:8000
python 01_basic_usage.py
python 02_structure_optimization.py

# 8. 监控（持续运行）
# 终端 1: 监控 GPU
watch -n 1 nvidia-smi

# 终端 2: 监控容器
docker stats

# 终端 3: 监控日志
./scripts/deploy.sh logs
```

**部署成功标志：**
- ✅ 所有 Worker 健康检查通过
- ✅ Gateway 可访问 http://localhost:8000/docs
- ✅ `/api/v1/models` 返回所有模型
- ✅ nvidia-smi 显示 5 个进程各占用一张 GPU
- ✅ examples 可以成功运行
