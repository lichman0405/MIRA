# MIRA Gateway API 文档

## 概述

MIRA Gateway 提供统一的 API 入口，自动路由请求到对应的 Worker 服务。

| 项目 | 信息 |
|------|------|
| **基础 URL** | `http://localhost:8000` |
| **API 前缀** | `/api/v1` |
| **数据格式** | JSON |
| **超时建议** | 10-30 分钟（计算任务） |

### 交互式文档

- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

---

## 健康检查

### GET /health

检查服务状态和 Worker 可用性。

**响应示例:**
```json
{
  "status": "healthy",
  "gateway": "healthy",
  "workers": {
    "mace-orb": {
      "status": "healthy",
      "service": "MACE-ORB",
      "models": ["mace-mp", "mace-off23", "orb-v2", "orb-v3"],
      "device": "cuda",
      "gpu_name": "NVIDIA RTX 4090"
    },
    "fairchem-sevennet": {
      "status": "unavailable",
      "error": "Connection refused"
    }
  },
  "available_workers": 1,
  "total_workers": 5
}
```

---

## 模型接口

### GET /api/v1/models

列出所有可用模型。

**响应示例:**
```json
{
  "total": 8,
  "models": [
    {"name": "mace-mp", "service": "mace-orb", "url": "http://mace-orb:8001"},
    {"name": "orb-v2", "service": "mace-orb", "url": "http://mace-orb:8001"}
  ]
}
```

### GET /api/v1/services

列出所有服务状态。

---

## 计算接口

所有计算接口使用相同的原子数据格式：

```json
{
  "atoms": {
    "symbols": ["Cu", "Cu", "Cu", "Cu"],
    "positions": [[0,0,0], [1.8,1.8,0], [1.8,0,1.8], [0,1.8,1.8]],
    "cell": [[3.6,0,0], [0,3.6,0], [0,0,3.6]],
    "pbc": [true, true, true]
  },
  "model_name": "mace-mp"
}
```

### POST /api/v1/single_point

单点能量计算。

**请求体:**
```json
{
  "atoms": { ... },
  "model_name": "mace-mp",
  "compute_stress": true,
  "compute_forces": true
}
```

**响应示例:**
```json
{
  "energy": -14.256,
  "energy_per_atom": -3.564,
  "forces": [[0.01, 0.02, -0.01], ...],
  "stress": [0.1, 0.1, 0.1, 0.0, 0.0, 0.0],
  "max_force": 0.025
}
```

### POST /api/v1/optimization

结构优化。

**请求体:**
```json
{
  "atoms": { ... },
  "model_name": "mace-mp",
  "fmax": 0.05,
  "max_steps": 500,
  "optimizer": "BFGS",
  "use_d3": true,
  "fix_cell": false
}
```

**参数说明:**
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| fmax | float | 0.05 | 力收敛标准 (eV/Å) |
| max_steps | int | 500 | 最大优化步数 |
| optimizer | str | "BFGS" | 优化器: BFGS/FIRE/LBFGS |
| use_d3 | bool | true | 启用 D3 色散校正 |
| fix_cell | bool | false | 固定晶胞 |

**响应示例:**
```json
{
  "initial_energy": -14.256,
  "final_energy": -14.892,
  "energy_change": -0.636,
  "steps": 42,
  "converged": true,
  "final_atoms": { ... }
}
```

### POST /api/v1/stability

MD 稳定性测试。

**请求体:**
```json
{
  "atoms": { ... },
  "model_name": "mace-mp",
  "temperature": 300.0,
  "pressure": 0.0,
  "timestep": 1.0,
  "equilibration_steps": 1000,
  "production_steps": 5000,
  "use_d3": true
}
```

**参数说明:**
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| temperature | float | 300.0 | 温度 (K) |
| pressure | float | 0.0 | 压力 (GPa), 0=NVT |
| timestep | float | 1.0 | 时间步长 (fs) |
| equilibration_steps | int | 1000 | 平衡步数 |
| production_steps | int | 5000 | 生产步数 |

### POST /api/v1/bulk_modulus

体积模量计算。

**请求体:**
```json
{
  "atoms": { ... },
  "model_name": "mace-mp",
  "strain_range": 0.06,
  "num_points": 7,
  "use_d3": true
}
```

**响应示例:**
```json
{
  "bulk_modulus": 142.5,
  "equilibrium_volume": 47.2,
  "bulk_modulus_derivative": 4.2,
  "r_squared": 0.9998
}
```

### POST /api/v1/heat_capacity

热容计算（声子）。

**请求体:**
```json
{
  "atoms": { ... },
  "model_name": "mace-mp",
  "temperatures": [100, 200, 300, 400, 500],
  "supercell": [2, 2, 2],
  "use_d3": true
}
```

### POST /api/v1/multi_model

多模型并行计算。

**请求体:**
```json
{
  "atoms": { ... },
  "model_names": ["mace-mp", "orb-v2", "m3gnet"],
  "task": "single_point",
  "task_params": {}
}
```

---

## Python 客户端示例

```python
import requests
from ase.io import read

# 配置
BASE_URL = "http://localhost:8000/api/v1"

# 加载结构
atoms = read("structure.cif")

# 转换为 API 格式
atoms_data = {
    "symbols": list(atoms.get_chemical_symbols()),
    "positions": atoms.get_positions().tolist(),
    "cell": atoms.get_cell().tolist(),
    "pbc": [bool(p) for p in atoms.get_pbc()]
}

# 单点能量计算
response = requests.post(
    f"{BASE_URL}/single_point",
    json={
        "atoms": atoms_data,
        "model_name": "mace-mp"
    },
    timeout=300
)
print(response.json())

# 结构优化
response = requests.post(
    f"{BASE_URL}/optimization",
    json={
        "atoms": atoms_data,
        "model_name": "mace-mp",
        "fmax": 0.05,
        "max_steps": 200
    },
    timeout=600
)
print(response.json())
```

---

## 错误处理

| HTTP 状态码 | 说明 |
|-------------|------|
| 200 | 成功 |
| 400 | 请求参数错误 |
| 404 | 模型未找到 |
| 503 | Worker 服务不可用 |
| 504 | 计算超时 |

**错误响应格式:**
```json
{
  "detail": "Model 'unknown-model' not found"
}
```
