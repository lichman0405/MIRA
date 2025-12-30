# GRACE Worker Service
# 端口: 8004
# 注意: GRACE 使用 TensorFlow，与 PyTorch 模型隔离

FROM nvidia/cuda:12.1-cudnn8-runtime-ubuntu22.04

# 环境变量
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    DEBIAN_FRONTEND=noninteractive \
    TF_CPP_MIN_LOG_LEVEL=2

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
    python3.10 \
    python3.10-dev \
    python3.10-venv \
    python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1 \
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1

# 升级 pip
RUN python -m pip install --upgrade pip setuptools wheel

# 安装 TensorFlow (GRACE 依赖)
RUN pip install tensorflow[and-cuda]

# 安装公共依赖
RUN pip install \
    numpy>=1.24.0 \
    scipy>=1.10.0 \
    ase>=3.27.0 \
    phonopy>=2.20.0 \
    fastapi>=0.104.0 \
    uvicorn[standard]>=0.24.0 \
    pydantic>=2.5.0 \
    httpx>=0.25.0

# 安装 GRACE (tensorpotential)
RUN pip install tensorpotential

# 安装 PyTorch CPU 版本 (仅用于 torch_dftd)
RUN pip install torch --index-url https://download.pytorch.org/whl/cpu
RUN pip install torch-dftd>=0.4.0

# 创建工作目录
WORKDIR /app

# 复制共享模块
COPY services/shared /app/shared

# 复制服务代码
COPY services/grace /app/grace_service

# 暴露端口
EXPOSE 8004

# 健康检查
HEALTHCHECK --interval=30s --timeout=30s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8004/health || exit 1

# 启动命令
CMD ["python", "-m", "uvicorn", "grace_service.main:app", "--host", "0.0.0.0", "--port", "8004"]
