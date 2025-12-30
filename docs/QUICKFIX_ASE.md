# ASE 导入错误快速修复指南

## 错误信息
```
ImportError: cannot import name 'ExpCellFilter' from 'ase.constraints'
```

## 原因
ASE >= 3.23 版本中，`ExpCellFilter` 和 `StrainFilter` 从 `ase.constraints` 移到了 `ase.filters`

## 已修复文件
✅ `services/shared/worker_base.py`
✅ `app/services/optimization.py`
✅ 所有 Worker Dockerfile 已更新注释

## 立即修复步骤

### 在测试服务器上（192.168.100.207）

```bash
# 1. 拉取最新代码
cd MIRA
git pull

# 2. 重新构建 CPU 镜像
./scripts/deploy.sh build-cpu

# 3. 停止旧容器
./scripts/deploy.sh stop

# 4. 启动新容器
./scripts/deploy.sh test-cpu

# 5. 检查日志（应该没有 ImportError 了）
docker-compose -f deployments/docker-compose.cpu-test.yml logs mace-orb-worker | head -50

# 6. 测试 API
curl http://localhost:8000/health
curl http://localhost:8000/api/v1/models
```

### 如果还有问题

```bash
# 进入容器检查
docker exec -it mira-mace-orb-worker bash

# 在容器内测试导入
python -c "from ase.filters import ExpCellFilter; print('✓ 导入成功')"

# 检查 ASE 版本
python -c "import ase; print(f'ASE 版本: {ase.__version__}')"

# 运行完整测试
python /app/scripts/test_ase_imports.py
```

## 验证修复成功

访问 Gateway 查看可用模型：
```bash
curl http://192.168.100.207:8000/api/v1/models
```

应该看到：
```json
{
  "available_models": [
    "mace-mp",
    "orb-v2",
    ...
  ],
  "model_count": 5
}
```

## 相关文档
- 详细说明：[docs/ASE_COMPATIBILITY_FIX.md](ASE_COMPATIBILITY_FIX.md)
- API 文档：[docs/API.md](API.md)
- CPU 兼容性：[docs/CPU_COMPATIBILITY.md](CPU_COMPATIBILITY.md)
