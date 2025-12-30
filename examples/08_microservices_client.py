#!/usr/bin/env python3
"""
MIRA 微服务客户端示例

演示如何调用 MIRA Gateway 进行计算

运行前确保:
1. 服务已启动: 
   - 本地: uvicorn app.main:app --host 0.0.0.0 --port 8000
   - Docker: ./scripts/deploy.sh test

支持环境变量:
    MIRA_GATEWAY_URL=http://192.168.100.207:8000  # 测试服务器
"""
import os
import httpx
import asyncio
from typing import List, Dict, Any

# ========== 配置 ==========
GATEWAY_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000")


class MIRAClient:
    """MIRA 微服务客户端"""
    
    def __init__(self, base_url: str = None):
        self.base_url = base_url or GATEWAY_URL
        self.client = httpx.AsyncClient(base_url=self.base_url, timeout=600.0)
        print(f"[配置] MIRA Gateway: {self.base_url}")
    
    async def close(self):
        await self.client.aclose()
    
    async def health(self) -> Dict[str, Any]:
        """健康检查"""
        response = await self.client.get("/health")
        response.raise_for_status()
        return response.json()
    
    async def list_models(self) -> List[Dict[str, Any]]:
        """列出可用模型"""
        response = await self.client.get("/api/v1/models")
        response.raise_for_status()
        return response.json()["models"]
    
    async def single_point(
        self,
        model_name: str,
        symbols: List[str],
        positions: List[List[float]],
        cell: List[List[float]] = None,
        pbc: List[bool] = None
    ) -> Dict[str, Any]:
        """单点能量计算"""
        payload = {
            "model_name": model_name,
            "atoms": {
                "symbols": symbols,
                "positions": positions,
                "cell": cell,
                "pbc": pbc or [True, True, True]
            }
        }
        response = await self.client.post("/api/v1/single_point", json=payload)
        response.raise_for_status()
        return response.json()
    
    async def optimization(
        self,
        model_name: str,
        symbols: List[str],
        positions: List[List[float]],
        cell: List[List[float]] = None,
        fmax: float = 0.05,
        max_steps: int = 500,
        use_d3: bool = True
    ) -> Dict[str, Any]:
        """结构优化"""
        payload = {
            "model_name": model_name,
            "atoms": {
                "symbols": symbols,
                "positions": positions,
                "cell": cell,
                "pbc": [True, True, True]
            },
            "fmax": fmax,
            "max_steps": max_steps,
            "use_d3": use_d3
        }
        response = await self.client.post("/api/v1/optimization", json=payload)
        response.raise_for_status()
        return response.json()
    
    async def multi_model_compare(
        self,
        model_names: List[str],
        task: str,
        atoms: Dict[str, Any],
        task_params: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """多模型比较"""
        payload = {
            "model_names": model_names,
            "task": task,
            "atoms": atoms,
            "task_params": task_params or {}
        }
        response = await self.client.post("/api/v1/multi_model", json=payload)
        response.raise_for_status()
        return response.json()


async def demo():
    """演示 MIRA 客户端使用"""
    
    print("=" * 60)
    print("MIRA 微服务客户端演示")
    print("=" * 60)
    
    client = MIRAClient()
    
    try:
        # 1. 健康检查
        print("\n1. 健康检查...")
        health = await client.health()
        print(f"   Gateway: {health['gateway']}")
        print(f"   可用 Workers: {health['available_workers']}/{health['total_workers']}")
        
        # 2. 列出可用模型
        print("\n2. 可用模型...")
        models = await client.list_models()
        for model in models[:5]:  # 只显示前 5 个
            print(f"   - {model['name']} ({model['service']})")
        if len(models) > 5:
            print(f"   ... 共 {len(models)} 个模型")
        
        # 3. 单点能量计算
        print("\n3. 单点能量计算 (Cu FCC)...")
        
        # Cu FCC 结构
        cu_symbols = ["Cu", "Cu", "Cu", "Cu"]
        cu_positions = [
            [0.0, 0.0, 0.0],
            [1.8, 1.8, 0.0],
            [1.8, 0.0, 1.8],
            [0.0, 1.8, 1.8]
        ]
        cu_cell = [
            [3.6, 0.0, 0.0],
            [0.0, 3.6, 0.0],
            [0.0, 0.0, 3.6]
        ]
        
        if models:
            model_name = models[0]["name"]
            print(f"   使用模型: {model_name}")
            
            result = await client.single_point(
                model_name=model_name,
                symbols=cu_symbols,
                positions=cu_positions,
                cell=cu_cell
            )
            print(f"   能量: {result['energy']:.4f} eV")
            print(f"   每原子能量: {result['energy_per_atom']:.4f} eV/atom")
        
        # 4. 多模型比较
        print("\n4. 多模型比较...")
        
        # 找出可用的模型
        available_model_names = [m["name"] for m in models[:3]]  # 取前 3 个
        
        if len(available_model_names) >= 2:
            print(f"   比较模型: {available_model_names}")
            
            comparison = await client.multi_model_compare(
                model_names=available_model_names,
                task="single_point",
                atoms={
                    "symbols": cu_symbols,
                    "positions": cu_positions,
                    "cell": cu_cell,
                    "pbc": [True, True, True]
                }
            )
            
            print("\n   比较结果:")
            for r in comparison["results"]:
                if r["status"] == "success":
                    print(f"   - {r['model']}: {r['result']['energy']:.4f} eV")
                else:
                    print(f"   - {r['model']}: 错误 - {r['error']}")
        
        print("\n" + "=" * 60)
        print("演示完成!")
        print("=" * 60)
        
    except httpx.ConnectError:
        print("错误: 无法连接到 MIRA Gateway")
        print("请确保服务已启动: ./scripts/deploy.sh test")
    except Exception as e:
        print(f"错误: {e}")
    finally:
        await client.close()


if __name__ == "__main__":
    asyncio.run(demo())
