# MIRA - MOF ML Force Field Benchmark Service
# Dockerfile for GPU-enabled container

# 基础镜像：NVIDIA CUDA with cuDNN
FROM nvidia/cuda:12.1-cudnn8-devel-ubuntu22.04

# 环境变量
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# 安装系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# 设置 Python 3.11 为默认
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# 创建工作目录
WORKDIR /app

# 升级 pip
RUN python -m pip install --upgrade pip setuptools wheel

# 复制依赖文件
COPY requirements.txt .

# 安装 Python 依赖（分层安装以利用缓存）
# 先安装 PyTorch
RUN pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 安装 ASE 和基础科学计算包
RUN pip install ase>=3.27.0 numpy scipy phonopy

# 安装 FastAPI 和相关包
RUN pip install fastapi uvicorn[standard] pydantic pydantic-settings python-multipart

# 安装 ML 模型相关包
RUN pip install torch-dftd e3nn

# 安装各模型包（可能需要较长时间）
RUN pip install mace-torch || echo "MACE installation skipped"
RUN pip install orb-models || echo "ORB installation skipped"
RUN pip install fairchem-core || echo "FAIRChem installation skipped"
RUN pip install sevenn || echo "SevenNet installation skipped"
RUN pip install matgl || echo "MatGL installation skipped"
RUN pip install huggingface_hub || echo "huggingface_hub installation skipped"

# 复制应用代码
COPY . .

# 创建数据目录
RUN mkdir -p /data/structures /data/results /data/models

# 设置环境变量
ENV STRUCTURES_DIR=/data/structures \
    RESULTS_DIR=/data/results \
    MODELS_CACHE_DIR=/data/models \
    DEFAULT_DEVICE=cuda

# 暴露端口
EXPOSE 8000

# 健康检查
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/api/v1/health || exit 1

# 启动命令
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
