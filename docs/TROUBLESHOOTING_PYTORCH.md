# PyTorch/torchvision 版本不兼容问题修复

## 错误信息

```
RuntimeError: operator torchvision::nms does not exist
```

## 问题原因

PyTorch 和 torchvision 版本不匹配导致。这个问题通常发生在：
- 分别安装了不同版本的 PyTorch 和 torchvision
- pip 自动解析依赖时选择了不兼容的版本组合

## 快速修复（推荐）

### 方法 1: 使用自动修复脚本

```bash
cd MIRA
python scripts/fix_torch_versions.py
```

脚本会：
1. 检查当前版本
2. 卸载旧版本（torch, torchvision, torchaudio, torchmetrics, mace-torch）
3. 安装兼容的新版本
4. 验证安装

### 方法 2: 手动修复

#### 步骤 1: 卸载旧版本

```bash
pip uninstall -y mace-torch torchmetrics torchvision torchaudio torch
```

#### 步骤 2: 安装兼容版本

**选项 A: PyTorch 2.5.1 + CUDA 12.4 (推荐)**
```bash
pip install torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 \
    --index-url https://download.pytorch.org/whl/cu124
```

**选项 B: PyTorch 2.4.1 + CUDA 12.1 (稳定)**
```bash
pip install torch==2.4.1 torchvision==0.19.1 torchaudio==2.4.1 \
    --index-url https://download.pytorch.org/whl/cu121
```

**选项 C: CPU 版本（无 GPU）**
```bash
pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cpu
```

#### 步骤 3: 重新安装依赖

```bash
pip install torchmetrics>=1.0.0
pip install mace-torch
```

#### 步骤 4: 验证安装

```bash
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import torchvision; print(f'torchvision: {torchvision.__version__}')"
python -c "from mace.calculators import mace_mp; print('✓ MACE 导入成功')"
```

## 兼容版本对照表

| PyTorch | torchvision | torchaudio | CUDA | 状态 |
|---------|-------------|------------|------|------|
| 2.5.1 | 0.20.1 | 2.5.1 | 12.4 | ✅ 推荐 |
| 2.4.1 | 0.19.1 | 2.4.1 | 12.1 | ✅ 稳定 |
| 2.3.1 | 0.18.1 | 2.3.1 | 12.1 | ✅ 可用 |
| 2.2.2 | 0.17.2 | 2.2.2 | 12.1 | ✅ 可用 |

## 测试修复

运行 examples 验证：

```bash
cd MIRA
python examples/01_basic_usage.py
```

应该看到：
```
cuequivariance or cuequivariance_torch is not available. Cuequivariance acceleration will be disabled.
✓ 检查通过
✓ 连接到 MIRA 网关
...
```

## 为什么 Docker 容器没问题？

Docker 镜像在构建时固定了兼容的版本：

```dockerfile
# docker/Dockerfile.mace-orb-cpu (举例)
RUN pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cpu
```

Docker 构建过程会自动选择兼容的版本组合，而本地环境可能因为：
- 之前安装过不同版本的包
- pip 缓存问题
- 依赖解析顺序不同

导致版本不兼容。

## 进一步排查

如果修复后仍有问题，检查：

### 1. 查看所有 PyTorch 相关包版本

```bash
pip list | grep -E "torch|vision|audio|metrics"
```

### 2. 清理 pip 缓存

```bash
pip cache purge
```

### 3. 完全重建虚拟环境

```bash
# 备份当前环境（可选）
pip freeze > old_requirements.txt

# 删除旧环境
deactivate
rm -rf venv

# 创建新环境
python -m venv venv
source venv/bin/activate  # Linux/macOS
# 或: venv\Scripts\activate  # Windows

# 安装基础包
pip install --upgrade pip setuptools wheel

# 安装兼容的 PyTorch
pip install torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 \
    --index-url https://download.pytorch.org/whl/cu124

# 安装其他依赖
pip install -r requirements.txt
python scripts/install_models.py --combo-a
```

## cuequivariance 警告

```
cuequivariance or cuequivariance_torch is not available. Cuequivariance acceleration will be disabled.
```

这是**正常的警告**，不影响功能。cuequivariance 是 MACE 的可选加速库，不安装也能正常运行。

## 相关资源

- [PyTorch 官方安装指南](https://pytorch.org/get-started/locally/)
- [torchvision 兼容性矩阵](https://github.com/pytorch/vision#installation)
- [MACE 安装文档](https://github.com/ACEsuit/mace)

## 总结

问题根源是 **PyTorch 和 torchvision 版本不匹配**。

**推荐修复方案：**
1. 运行 `python scripts/fix_torch_versions.py`
2. 选择与你的 CUDA 版本匹配的选项
3. 重新测试 examples

**预防措施：**
- 安装 PyTorch 时使用官方推荐的命令
- 同时安装 torch, torchvision, torchaudio
- 使用虚拟环境隔离项目依赖
