"""
MIRA - MOF ML Force Field Benchmark Service

FastAPI 主入口
"""
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import get_settings
from .api.v1.router import api_router
from .dependencies import startup_handler, shutdown_handler


@asynccontextmanager
async def lifespan(app: FastAPI):
    """应用生命周期管理"""
    # 启动
    await startup_handler()
    yield
    # 关闭
    await shutdown_handler()


def create_app() -> FastAPI:
    """创建 FastAPI 应用"""
    settings = get_settings()
    
    app = FastAPI(
        title="MIRA - MiQroEra Interatomic-potential Reliability Arena",
        description="""
## MIRA: MiQroEra Interatomic-potential Reliability Arena

A comprehensive RESTful API for benchmarking machine learning interatomic potentials on Metal-Organic Frameworks (MOFs).

**Author:** Shibo Li (shadow.li981@gmail.com)

### Features

- **Multiple Model Families**: Support for 20+ ML force field models including:
  - MACE (MP, OFF23, OMAT, MPA, ANI)
  - ORB (v2, v3, OMAT-v3-LoRA)
  - OMAT24 / FAIRChem
  - GRACE-2L / GRACE-2M
  - MatterSim
  - SevenNet (0, MF-ompa, l3i5)
  - PosEGNN

- **Computational Tasks**:
  - Structure Optimization (BFGS/FIRE/LBFGS with FrechetCellFilter)
  - MD Stability Testing (NVT/NPT dynamics)
  - Bulk Modulus Calculation (Birch-Murnaghan EOS fitting)
  - Heat Capacity (Phonon calculations with Phonopy)
  - QMOF Energy Evaluation
  - Interaction Energy Analysis

### Usage

1. Upload a structure (CIF/POSCAR/XYZ format)
2. Select a model from the available options
3. Submit a task for computation
4. Monitor progress and retrieve results

### Notes

- This service requires GPU for optimal performance
- ASE >= 3.27.0 is required for NPT dynamics support
- Large MD simulations may take significant time
        """,
        version="1.0.0",
        docs_url="/docs",
        redoc_url="/redoc",
        lifespan=lifespan
    )
    
    # CORS 中间件
    if settings.cors_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=settings.cors_origins,
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"]
        )
    
    # 注册 API 路由
    app.include_router(api_router, prefix="/api")
    
    # 根路径
    @app.get("/")
    async def root():
        return {
            "service": "MIRA",
            "description": "MiQroEra Interatomic-potential Reliability Arena",
            "author": "Shibo Li",
            "email": "shadow.li981@gmail.com",
            "version": "1.0.0",
            "docs": "/docs",
            "api": "/api/v1"
        }
    
    return app


# 创建应用实例
app = create_app()


if __name__ == "__main__":
    import uvicorn
    
    settings = get_settings()
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.debug
    )
