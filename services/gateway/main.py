#!/usr/bin/env python3
"""
MIRA Gateway Service

统一入口，将请求路由到对应的模型 Worker 服务

架构:
  用户 -> Gateway (8000) -> Worker Services (8001-8005)
"""
import os
import asyncio
from typing import Dict, List, Optional, Any
from contextlib import asynccontextmanager

import httpx
from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

# ============================================
# 配置
# ============================================
class ServiceConfig:
    """服务配置"""
    
    # Worker 服务地址 (Docker 网络内部)
    SERVICES = {
        "mace-orb": os.getenv("MACE_ORB_URL", "http://mira-mace-orb:8001"),
        "fairchem-sevennet": os.getenv("FAIRCHEM_SEVENNET_URL", "http://mira-fairchem-sevennet:8002"),
        "matgl": os.getenv("MATGL_URL", "http://mira-matgl:8003"),
        "grace": os.getenv("GRACE_URL", "http://mira-grace:8004"),
        "mattersim": os.getenv("MATTERSIM_URL", "http://mira-mattersim:8005"),
    }
    
    # 模型到服务的映射
    MODEL_SERVICE_MAP = {
        # MACE 系列
        "mace-mp": "mace-orb",
        "mace-off23": "mace-orb",
        "mace-omat": "mace-orb",
        "mace-mpa": "mace-orb",
        "mace-ani": "mace-orb",
        # ORB 系列
        "orb-v2": "mace-orb",
        "orb-d3-v2": "mace-orb",
        "orb-v3": "mace-orb",
        # FAIRChem 系列
        "omat24-base": "fairchem-sevennet",
        "omat24-large": "fairchem-sevennet",
        "eqv2-omat": "fairchem-sevennet",
        "eqv2-mptrj": "fairchem-sevennet",
        # SevenNet 系列
        "sevennet-0": "fairchem-sevennet",
        "sevennet-mf-ompa": "fairchem-sevennet",
        "sevennet-l3i5": "fairchem-sevennet",
        # MatGL 系列
        "m3gnet": "matgl",
        "m3gnet-direct": "matgl",
        "chgnet": "matgl",
        "chgnet-mpp": "matgl",
        # GRACE 系列
        "grace-2l": "grace",
        "grace-2l-omat": "grace",
        "grace-2m": "grace",
        # MatterSim
        "mattersim-5m": "mattersim",
    }
    
    @classmethod
    def get_service_url(cls, model_name: str) -> str:
        """获取模型对应的服务 URL"""
        service_key = cls.MODEL_SERVICE_MAP.get(model_name)
        if not service_key:
            raise ValueError(f"Unknown model: {model_name}")
        return cls.SERVICES[service_key]


# ============================================
# HTTP 客户端
# ============================================
class ServiceClient:
    """Worker 服务客户端"""
    
    def __init__(self):
        self.client = httpx.AsyncClient(timeout=600.0)  # 10 分钟超时
        self._service_status: Dict[str, bool] = {}
    
    async def close(self):
        await self.client.aclose()
    
    async def check_service(self, service_name: str, url: str) -> Dict[str, Any]:
        """检查服务健康状态"""
        try:
            response = await self.client.get(f"{url}/health", timeout=10.0)
            if response.status_code == 200:
                self._service_status[service_name] = True
                return response.json()
            else:
                self._service_status[service_name] = False
                return {"status": "unhealthy", "error": response.text}
        except Exception as e:
            self._service_status[service_name] = False
            return {"status": "unavailable", "error": str(e)}
    
    async def check_all_services(self) -> Dict[str, Any]:
        """检查所有服务状态"""
        results = {}
        for name, url in ServiceConfig.SERVICES.items():
            results[name] = await self.check_service(name, url)
        return results
    
    async def get_available_models(self) -> List[Dict[str, Any]]:
        """获取所有可用模型"""
        models = []
        for name, url in ServiceConfig.SERVICES.items():
            try:
                response = await self.client.get(f"{url}/models", timeout=10.0)
                if response.status_code == 200:
                    data = response.json()
                    for model in data.get("models", []):
                        models.append({
                            "name": model,
                            "service": name,
                            "url": url
                        })
            except Exception:
                pass
        return models
    
    async def forward_request(
        self,
        model_name: str,
        endpoint: str,
        payload: Dict[str, Any]
    ) -> Dict[str, Any]:
        """转发请求到对应的 Worker"""
        url = ServiceConfig.get_service_url(model_name)
        
        try:
            response = await self.client.post(
                f"{url}/{endpoint}",
                json=payload,
                timeout=600.0
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                raise HTTPException(
                    status_code=response.status_code,
                    detail=response.json().get("detail", response.text)
                )
        except httpx.TimeoutException:
            raise HTTPException(status_code=504, detail="Worker service timeout")
        except httpx.ConnectError:
            raise HTTPException(status_code=503, detail=f"Worker service unavailable")


# ============================================
# 请求/响应模型
# ============================================
class AtomData(BaseModel):
    symbols: List[str]
    positions: List[List[float]]
    cell: Optional[List[List[float]]] = None
    pbc: List[bool] = [True, True, True]


class SinglePointRequest(BaseModel):
    atoms: AtomData
    model_name: str
    compute_stress: bool = True
    compute_forces: bool = True


class OptimizationRequest(BaseModel):
    atoms: AtomData
    model_name: str
    fmax: float = 0.05
    max_steps: int = 500
    optimizer: str = "BFGS"
    use_d3: bool = True
    fix_cell: bool = False


class StabilityRequest(BaseModel):
    atoms: AtomData
    model_name: str
    temperature: float = 300.0
    pressure: float = 0.0
    timestep: float = 1.0
    equilibration_steps: int = 1000
    production_steps: int = 5000
    use_d3: bool = True


class BulkModulusRequest(BaseModel):
    atoms: AtomData
    model_name: str
    strain_range: float = 0.06
    num_points: int = 7
    use_d3: bool = True


class HeatCapacityRequest(BaseModel):
    atoms: AtomData
    model_name: str
    temperatures: List[float] = [100, 200, 300, 400, 500]
    supercell: List[int] = [2, 2, 2]
    use_d3: bool = True


class MultiModelRequest(BaseModel):
    """多模型比较请求"""
    atoms: AtomData
    model_names: List[str]
    task: str = "single_point"
    task_params: Dict[str, Any] = {}


# ============================================
# FastAPI 应用
# ============================================
client = ServiceClient()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """应用生命周期"""
    print("MIRA Gateway starting...")
    print(f"Configured services: {list(ServiceConfig.SERVICES.keys())}")
    yield
    await client.close()
    print("MIRA Gateway stopped")


app = FastAPI(
    title="MIRA - MOF ML Force Field Benchmark Service",
    description="""
MiQroEra Interatomic-potential Reliability Arena

Gateway 服务 - 统一入口，自动路由到各模型 Worker 服务

## 架构

```
用户请求 -> Gateway (8000) -> Worker Services
                              ├── MACE+ORB (8001)
                              ├── FAIRChem+SevenNet (8002)
                              ├── MatGL (8003)
                              ├── GRACE (8004)
                              └── MatterSim (8005)
```

## 支持的模型

- **MACE**: mace-mp, mace-off23, mace-omat, mace-mpa, mace-ani
- **ORB**: orb-v2, orb-d3-v2, orb-v3
- **FAIRChem**: omat24-base, omat24-large, eqv2-omat, eqv2-mptrj
- **SevenNet**: sevennet-0, sevennet-mf-ompa, sevennet-l3i5
- **MatGL**: m3gnet, chgnet
- **GRACE**: grace-2l, grace-2l-omat, grace-2m
- **MatterSim**: mattersim-5m
    """,
    version="2.0.0",
    lifespan=lifespan
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ============================================
# API 路由
# ============================================
@app.get("/")
async def root():
    """根路径"""
    return {
        "service": "MIRA Gateway",
        "version": "2.0.0",
        "description": "MiQroEra Interatomic-potential Reliability Arena"
    }


@app.get("/health")
async def health():
    """健康检查"""
    services = await client.check_all_services()
    healthy_count = sum(1 for s in services.values() if s.get("status") == "healthy")
    
    return {
        "status": "healthy" if healthy_count > 0 else "degraded",
        "gateway": "healthy",
        "workers": services,
        "available_workers": healthy_count,
        "total_workers": len(services)
    }


@app.get("/api/v1/models")
async def list_models():
    """列出所有可用模型"""
    models = await client.get_available_models()
    return {
        "total": len(models),
        "models": models
    }


@app.get("/api/v1/services")
async def list_services():
    """列出所有服务状态"""
    return await client.check_all_services()


@app.post("/api/v1/single_point")
async def single_point(request: SinglePointRequest):
    """单点能量计算"""
    return await client.forward_request(
        request.model_name,
        "single_point",
        request.model_dump()
    )


@app.post("/api/v1/optimization")
async def optimization(request: OptimizationRequest):
    """结构优化"""
    return await client.forward_request(
        request.model_name,
        "optimization",
        request.model_dump()
    )


@app.post("/api/v1/stability")
async def stability(request: StabilityRequest):
    """稳定性分析"""
    return await client.forward_request(
        request.model_name,
        "stability",
        request.model_dump()
    )


@app.post("/api/v1/bulk_modulus")
async def bulk_modulus(request: BulkModulusRequest):
    """体积模量计算"""
    return await client.forward_request(
        request.model_name,
        "bulk_modulus",
        request.model_dump()
    )


@app.post("/api/v1/heat_capacity")
async def heat_capacity(request: HeatCapacityRequest):
    """热容计算"""
    return await client.forward_request(
        request.model_name,
        "heat_capacity",
        request.model_dump()
    )


@app.post("/api/v1/multi_model")
async def multi_model_compare(request: MultiModelRequest):
    """多模型并行计算比较"""
    
    async def run_single_model(model_name: str) -> Dict[str, Any]:
        """对单个模型运行任务"""
        try:
            payload = {
                "atoms": request.atoms.model_dump(),
                "model_name": model_name,
                **request.task_params
            }
            result = await client.forward_request(
                model_name,
                request.task,
                payload
            )
            return {"model": model_name, "status": "success", "result": result}
        except HTTPException as e:
            return {"model": model_name, "status": "error", "error": e.detail}
        except Exception as e:
            return {"model": model_name, "status": "error", "error": str(e)}
    
    # 并行运行所有模型
    tasks = [run_single_model(model) for model in request.model_names]
    results = await asyncio.gather(*tasks)
    
    return {
        "task": request.task,
        "models_count": len(request.model_names),
        "results": results
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
