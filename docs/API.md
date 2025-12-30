# MIRA API 接口完整使用说明

## 目录

1. [概述](#概述)
2. [认证与基础配置](#认证与基础配置)
3. [模型管理接口](#模型管理接口)
4. [结构管理接口](#结构管理接口)
5. [任务管理接口](#任务管理接口)
6. [结果查询接口](#结果查询接口)
7. [完整使用流程示例](#完整使用流程示例)
8. [错误处理](#错误处理)
9. [最佳实践](#最佳实践)

---

## 概述

### 基础信息

| 项目 | 信息 |
|------|------|
| **基础 URL** | `http://localhost:8000` (本地) / `http://your-server:8000` (远程) |
| **API 版本** | v1 |
| **API 前缀** | `/api/v1` |
| **协议** | HTTP/HTTPS |
| **数据格式** | JSON |
| **编码** | UTF-8 |

### 交互式文档

MIRA 提供自动生成的交互式 API 文档：

- **Swagger UI**: `http://localhost:8000/docs`
- **ReDoc**: `http://localhost:8000/redoc`

### 健康检查

```bash
# 检查服务状态
curl http://localhost:8000/health
```

**响应示例：**
```json
{
  "status": "healthy",
  "version": "1.0.0",
  "timestamp": "2025-12-30T10:30:00Z"
}
```

---

## 认证与基础配置

### 环境变量配置

```bash
# 设置 Gateway 地址（客户端使用）
export MIRA_GATEWAY_URL=http://192.168.100.207:8000

# Windows PowerShell
$env:MIRA_GATEWAY_URL = "http://192.168.100.207:8000"
```

### Python 客户端配置

```python
import os
import requests

# 自动读取环境变量
GATEWAY_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000")
BASE_URL = f"{GATEWAY_URL}/api/v1"

# 设置超时时间
TIMEOUT = 600  # 10 分钟，适用于长时间计算任务
```

---

## 模型管理接口

### 1. 列出所有可用模型

**端点:** `GET /api/v1/models`

**描述:** 获取所有可用的 ML 力场模型列表

**查询参数:**
| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `family` | string | 否 | 按模型家族筛选 (mace/orb/matgl/grace等) |

**请求示例:**
```bash
# 获取所有模型
curl http://localhost:8000/api/v1/models

# 筛选 MACE 家族模型
curl http://localhost:8000/api/v1/models?family=mace
```

**Python 示例:**
```python
import requests

response = requests.get(f"{BASE_URL}/models")
models = response.json()

for model in models:
    print(f"{model['key']}: {model['description']}")
```

**响应示例:**
```json
[
  {
    "key": "mace-mp",
    "name": "MACE-MP-0",
    "family": "mace",
    "description": "Materials Project 2023 基础模型",
    "elements": ["H", "C", "N", "O", "F", "...", "Bi"],
    "supported_tasks": ["energy", "forces", "stress", "optimization", "md"],
    "references": "https://doi.org/10.48550/arXiv.2401.00096"
  },
  {
    "key": "orb-v2",
    "name": "ORB v2",
    "family": "orb",
    "description": "Orbital Materials v2 通用模型",
    "elements": ["all"],
    "supported_tasks": ["energy", "forces", "optimization", "md"],
    "references": "https://doi.org/..."
  }
]
```

---

### 2. 列出模型家族

**端点:** `GET /api/v1/models/families`

**描述:** 获取所有支持的模型家族列表

**请求示例:**
```bash
curl http://localhost:8000/api/v1/models/families
```

**响应示例:**
```json
["mace", "orb", "omat24", "grace", "mattersim", "sevennet", "matgl"]
```

---

### 3. 获取特定模型信息

**端点:** `GET /api/v1/models/{model_key}`

**描述:** 获取指定模型的详细信息

**路径参数:**
| 参数 | 类型 | 说明 |
|------|------|------|
| `model_key` | string | 模型标识符，如 `mace-mp`, `orb-v2` |

**请求示例:**
```bash
curl http://localhost:8000/api/v1/models/mace-mp
```

**Python 示例:**
```python
response = requests.get(f"{BASE_URL}/models/mace-mp")
model_info = response.json()
print(f"支持的元素: {', '.join(model_info['elements'])}")
```

---

### 4. 检查模型加载状态

**端点:** `GET /api/v1/models/{model_key}/status`

**描述:** 检查模型是否已加载到内存

**请求示例:**
```bash
curl http://localhost:8000/api/v1/models/mace-mp/status
```

**响应示例:**
```json
{
  "model_key": "mace-mp",
  "loaded": true,
  "device": "cuda:0",
  "memory_usage_mb": 2048.5
}
```

---

## 结构管理接口

### 1. 上传结构文件（表单方式）

**端点:** `POST /api/v1/structures/upload`

**描述:** 通过表单提交结构文件内容

**请求体 (Form Data):**
| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `name` | string | 是 | 结构名称 |
| `format` | string | 是 | 文件格式: `cif`, `poscar`, `xyz` |
| `content` | string | 是 | 文件内容（文本） |

**请求示例 (cURL):**
```bash
curl -X POST http://localhost:8000/api/v1/structures/upload \
  -F "name=MOF-5" \
  -F "format=cif" \
  -F "content=@structure.cif"
```

**Python 示例:**
```python
with open("MOF-5.cif", "r") as f:
    content = f.read()

response = requests.post(
    f"{BASE_URL}/structures/upload",
    data={
        "name": "MOF-5",
        "format": "cif",
        "content": content
    }
)

structure_info = response.json()
structure_id = structure_info["id"]
print(f"结构已上传，ID: {structure_id}")
```

**响应示例:**
```json
{
  "id": "550e8400-e29b-41d4-a716-446655440000",
  "name": "MOF-5",
  "format": "cif",
  "formula": "Zn4O13C24H12",
  "num_atoms": 53,
  "volume": 2876.42,
  "density": 0.593,
  "upload_time": "2025-12-30T10:30:00Z",
  "metadata": {
    "space_group": "Fm-3m",
    "a": 25.832,
    "b": 25.832,
    "c": 25.832,
    "alpha": 90.0,
    "beta": 90.0,
    "gamma": 90.0
  }
}
```

---

### 2. 上传结构文件（文件上传）

**端点:** `POST /api/v1/structures/upload/file`

**描述:** 直接上传结构文件

**请求体 (Multipart Form):**
| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `file` | file | 是 | 结构文件 |
| `name` | string | 否 | 自定义名称（默认使用文件名） |

**请求示例:**
```bash
curl -X POST http://localhost:8000/api/v1/structures/upload/file \
  -F "file=@MOF-5.cif" \
  -F "name=MOF-5-optimized"
```

**Python 示例:**
```python
with open("MOF-5.cif", "rb") as f:
    files = {"file": f}
    data = {"name": "MOF-5"}
    response = requests.post(
        f"{BASE_URL}/structures/upload/file",
        files=files,
        data=data
    )
structure_id = response.json()["id"]
```

---

### 3. 列出所有结构

**端点:** `GET /api/v1/structures`

**描述:** 获取已上传的结构列表

**查询参数:**
| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `skip` | integer | 否 | 跳过前 N 个结果（分页） |
| `limit` | integer | 否 | 返回最多 N 个结果（默认 100） |

**请求示例:**
```bash
# 获取所有结构
curl http://localhost:8000/api/v1/structures

# 分页查询
curl http://localhost:8000/api/v1/structures?skip=0&limit=10
```

**响应示例:**
```json
{
  "total": 5,
  "structures": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "MOF-5",
      "formula": "Zn4O13C24H12",
      "num_atoms": 53,
      "upload_time": "2025-12-30T10:30:00Z"
    }
  ]
}
```

---

### 4. 获取结构详细信息

**端点:** `GET /api/v1/structures/{structure_id}`

**描述:** 获取指定结构的完整信息

**路径参数:**
| 参数 | 类型 | 说明 |
|------|------|------|
| `structure_id` | string | 结构 UUID |

**请求示例:**
```bash
curl http://localhost:8000/api/v1/structures/550e8400-e29b-41d4-a716-446655440000
```

**Python 示例:**
```python
response = requests.get(f"{BASE_URL}/structures/{structure_id}")
structure = response.json()
print(f"化学式: {structure['formula']}")
print(f"原子数: {structure['num_atoms']}")
```

---

### 5. 预览结构

**端点:** `GET /api/v1/structures/{structure_id}/preview`

**描述:** 获取结构的原子坐标信息（用于可视化）

**请求示例:**
```bash
curl http://localhost:8000/api/v1/structures/550e8400-e29b-41d4-a716-446655440000/preview
```

**响应示例:**
```json
{
  "symbols": ["Zn", "Zn", "O", "O", "C", ...],
  "positions": [
    [0.0, 0.0, 0.0],
    [12.916, 12.916, 0.0],
    ...
  ],
  "cell": [
    [25.832, 0.0, 0.0],
    [0.0, 25.832, 0.0],
    [0.0, 0.0, 25.832]
  ],
  "pbc": [true, true, true]
}
```

---

### 6. 删除结构

**端点:** `DELETE /api/v1/structures/{structure_id}`

**描述:** 删除指定结构及相关任务结果

**请求示例:**
```bash
curl -X DELETE http://localhost:8000/api/v1/structures/550e8400-e29b-41d4-a716-446655440000
```

**响应示例:**
```json
{
  "message": "Structure deleted successfully",
  "deleted_tasks": 3
}
```

---

## 任务管理接口

### 1. 提交结构优化任务

**端点:** `POST /api/v1/tasks/optimization`

**描述:** 使用指定 ML 模型优化结构

**请求体 (JSON):**
```json
{
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "model_key": "mace-mp",
  "fmax": 0.01,
  "max_steps": 500,
  "optimizer": "BFGS",
  "use_filter": true,
  "enable_d3": true
}
```

**参数说明:**
| 字段 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `structure_id` | string | 是 | - | 结构 UUID |
| `model_key` | string | 是 | - | 模型标识符 |
| `fmax` | float | 否 | 0.05 | 力收敛标准 (eV/Å) |
| `max_steps` | integer | 否 | 500 | 最大优化步数 |
| `optimizer` | string | 否 | BFGS | 优化器: BFGS/FIRE/LBFGS |
| `use_filter` | boolean | 否 | false | 是否优化晶胞 |
| `enable_d3` | boolean | 否 | false | 启用 DFT-D3 色散校正 |

**Python 示例:**
```python
response = requests.post(
    f"{BASE_URL}/tasks/optimization",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "fmax": 0.01,
        "max_steps": 500,
        "optimizer": "BFGS",
        "use_filter": True,
        "enable_d3": True
    }
)

task = response.json()
task_id = task["task_id"]
print(f"任务已提交，ID: {task_id}")
```

**响应示例:**
```json
{
  "task_id": "7c9e6679-7425-40de-944b-e07fc1f90ae7",
  "task_type": "optimization",
  "status": "pending",
  "model_key": "mace-mp",
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "created_at": "2025-12-30T10:35:00Z",
  "estimated_time_minutes": 5
}
```

---

### 2. 提交 MD 稳定性测试任务

**端点:** `POST /api/v1/tasks/stability`

**描述:** 对结构进行分子动力学模拟，评估热稳定性

**请求体示例:**
```json
{
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "model_key": "mace-mp",
  "nvt_steps": 5000,
  "nvt_temperature": 300,
  "npt_steps": 10000,
  "npt_temperature": 300,
  "npt_pressure": 1.01325,
  "timestep": 1.0,
  "pre_optimize": true,
  "enable_d3": true
}
```

**参数说明:**
| 字段 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `structure_id` | string | 是 | - | 结构 UUID |
| `model_key` | string | 是 | - | 模型标识符 |
| `nvt_steps` | integer | 否 | 5000 | NVT 平衡步数 |
| `nvt_temperature` | float | 否 | 300 | NVT 温度 (K) |
| `npt_steps` | integer | 否 | 10000 | NPT 生产步数 |
| `npt_temperature` | float | 否 | 300 | NPT 温度 (K) |
| `npt_pressure` | float | 否 | 1.01325 | NPT 压力 (bar) |
| `timestep` | float | 否 | 1.0 | 时间步长 (fs) |
| `pre_optimize` | boolean | 否 | true | 模拟前先优化 |
| `enable_d3` | boolean | 否 | false | 启用 D3 校正 |

**Python 示例:**
```python
response = requests.post(
    f"{BASE_URL}/tasks/stability",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "nvt_steps": 5000,
        "npt_steps": 10000,
        "npt_temperature": 300
    },
    timeout=TIMEOUT
)
task_id = response.json()["task_id"]
```

---

### 3. 提交体积模量计算任务

**端点:** `POST /api/v1/tasks/bulk-modulus`

**描述:** 通过 E-V 曲线计算材料的体积模量

**请求体示例:**
```json
{
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "model_key": "mace-mp",
  "strain_range": 0.1,
  "num_points": 9,
  "pre_optimize": true,
  "enable_d3": true
}
```

**参数说明:**
| 字段 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `structure_id` | string | 是 | - | 结构 UUID |
| `model_key` | string | 是 | - | 模型标识符 |
| `strain_range` | float | 否 | 0.1 | 应变范围 (±10%) |
| `num_points` | integer | 否 | 9 | 采样点数 |
| `pre_optimize` | boolean | 否 | true | 计算前先优化 |
| `enable_d3` | boolean | 否 | false | 启用 D3 校正 |

---

### 4. 提交热容计算任务

**端点:** `POST /api/v1/tasks/heat-capacity`

**描述:** 使用 Phonopy 计算声子和热容

**请求体示例:**
```json
{
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "model_key": "mace-mp",
  "supercell": [2, 2, 2],
  "displacement": 0.01,
  "temperature_range": [0, 1000, 50],
  "pre_optimize": true
}
```

**参数说明:**
| 字段 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `structure_id` | string | 是 | - | 结构 UUID |
| `model_key` | string | 是 | - | 模型标识符 |
| `supercell` | array[int] | 否 | [2,2,2] | 超胞尺寸 |
| `displacement` | float | 否 | 0.01 | 原子位移 (Å) |
| `temperature_range` | array[float] | 否 | [0,1000,50] | [起始,终止,步长] K |
| `pre_optimize` | boolean | 否 | true | 计算前先优化 |

---

### 5. 查询任务列表

**端点:** `GET /api/v1/tasks`

**描述:** 获取所有任务列表

**查询参数:**
| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `status` | string | 否 | 筛选状态: pending/running/completed/failed |
| `task_type` | string | 否 | 筛选类型: optimization/stability/etc |
| `skip` | integer | 否 | 分页：跳过条数 |
| `limit` | integer | 否 | 分页：限制条数 |

**请求示例:**
```bash
# 查询所有正在运行的任务
curl http://localhost:8000/api/v1/tasks?status=running

# 查询优化任务
curl http://localhost:8000/api/v1/tasks?task_type=optimization&limit=10
```

**Python 示例:**
```python
response = requests.get(
    f"{BASE_URL}/tasks",
    params={"status": "completed", "limit": 10}
)
tasks = response.json()
```

---

### 6. 查询任务状态

**端点:** `GET /api/v1/tasks/{task_id}`

**描述:** 获取特定任务的详细状态

**请求示例:**
```bash
curl http://localhost:8000/api/v1/tasks/7c9e6679-7425-40de-944b-e07fc1f90ae7
```

**Python 示例（轮询直到完成）:**
```python
import time

while True:
    response = requests.get(f"{BASE_URL}/tasks/{task_id}")
    task = response.json()
    
    status = task["status"]
    print(f"任务状态: {status}")
    
    if status == "completed":
        print("任务完成！")
        break
    elif status == "failed":
        print(f"任务失败: {task.get('error_message')}")
        break
    
    time.sleep(5)  # 每 5 秒查询一次
```

**响应示例:**
```json
{
  "task_id": "7c9e6679-7425-40de-944b-e07fc1f90ae7",
  "task_type": "optimization",
  "status": "running",
  "model_key": "mace-mp",
  "structure_id": "550e8400-e29b-41d4-a716-446655440000",
  "created_at": "2025-12-30T10:35:00Z",
  "started_at": "2025-12-30T10:35:15Z",
  "progress": 42,
  "current_step": 210,
  "total_steps": 500,
  "estimated_remaining_minutes": 3
}
```

---

### 7. 取消任务

**端点:** `POST /api/v1/tasks/{task_id}/cancel`

**描述:** 取消正在运行或等待的任务

**请求示例:**
```bash
curl -X POST http://localhost:8000/api/v1/tasks/7c9e6679-7425-40de-944b-e07fc1f90ae7/cancel
```

**响应示例:**
```json
{
  "message": "Task cancelled successfully",
  "task_id": "7c9e6679-7425-40de-944b-e07fc1f90ae7"
}
```

---

## 结果查询接口

### 1. 获取任务结果

**端点:** `GET /api/v1/results/{task_id}`

**描述:** 获取已完成任务的详细结果

**请求示例:**
```bash
curl http://localhost:8000/api/v1/results/7c9e6679-7425-40de-944b-e07fc1f90ae7
```

**Python 示例:**
```python
response = requests.get(f"{BASE_URL}/results/{task_id}")
result = response.json()

if result["task_type"] == "optimization":
    print(f"最终能量: {result['final_energy']} eV")
    print(f"优化步数: {result['num_steps']}")
    print(f"收敛: {result['converged']}")
```

**优化任务响应示例:**
```json
{
  "task_id": "7c9e6679-7425-40de-944b-e07fc1f90ae7",
  "task_type": "optimization",
  "status": "completed",
  "model_key": "mace-mp",
  "initial_energy": -1250.3456,
  "final_energy": -1268.7890,
  "energy_change": -18.4434,
  "num_steps": 187,
  "converged": true,
  "final_fmax": 0.0087,
  "initial_volume": 2876.42,
  "final_volume": 2854.18,
  "volume_change_percent": -0.77,
  "rmsd": 0.234,
  "computation_time_seconds": 245.67,
  "optimized_structure": {
    "symbols": ["Zn", "Zn", ...],
    "positions": [[0.01, -0.02, 0.0], ...],
    "cell": [[25.745, 0.0, 0.0], ...]
  }
}
```

**MD 稳定性测试响应示例:**
```json
{
  "task_id": "...",
  "task_type": "stability",
  "status": "completed",
  "nvt_duration_ps": 5.0,
  "npt_duration_ps": 10.0,
  "avg_temperature": 299.8,
  "temp_std": 12.5,
  "final_rmsd": 0.456,
  "volume_drift_percent": 1.2,
  "stability_assessment": "stable",
  "trajectory_file": "/path/to/trajectory.traj"
}
```

---

### 2. 下载结果文件

**端点:** `GET /api/v1/results/{task_id}/download/{file_type}`

**描述:** 下载任务生成的文件

**路径参数:**
| 参数 | 类型 | 说明 |
|------|------|------|
| `task_id` | string | 任务 UUID |
| `file_type` | string | 文件类型: `structure`/`trajectory`/`log` |

**请求示例:**
```bash
# 下载优化后的结构
curl -O http://localhost:8000/api/v1/results/7c9e6679-7425-40de-944b-e07fc1f90ae7/download/structure

# 下载 MD 轨迹
curl -O http://localhost:8000/api/v1/results/7c9e6679-7425-40de-944b-e07fc1f90ae7/download/trajectory
```

**Python 示例:**
```python
# 下载优化后的结构
response = requests.get(
    f"{BASE_URL}/results/{task_id}/download/structure",
    stream=True
)

with open("optimized.cif", "wb") as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)
```

---

### 3. 获取结果摘要

**端点:** `GET /api/v1/results/{task_id}/summary`

**描述:** 获取结果的简要摘要（不包含完整数据）

**请求示例:**
```bash
curl http://localhost:8000/api/v1/results/7c9e6679-7425-40de-944b-e07fc1f90ae7/summary
```

**响应示例:**
```json
{
  "task_id": "7c9e6679-7425-40de-944b-e07fc1f90ae7",
  "task_type": "optimization",
  "status": "completed",
  "model_key": "mace-mp",
  "key_metrics": {
    "final_energy": -1268.7890,
    "converged": true,
    "num_steps": 187
  },
  "computation_time_seconds": 245.67,
  "completed_at": "2025-12-30T10:39:15Z"
}
```

---

## 完整使用流程示例

### 场景 1: 结构优化

```python
import requests
import time
import os

# 配置
BASE_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000") + "/api/v1"

# 1. 上传结构
with open("MOF-5.cif", "r") as f:
    content = f.read()

response = requests.post(
    f"{BASE_URL}/structures/upload",
    data={
        "name": "MOF-5",
        "format": "cif",
        "content": content
    }
)
structure_id = response.json()["id"]
print(f"✓ 结构已上传: {structure_id}")

# 2. 提交优化任务
response = requests.post(
    f"{BASE_URL}/tasks/optimization",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "fmax": 0.01,
        "max_steps": 500,
        "optimizer": "BFGS",
        "use_filter": True,
        "enable_d3": True
    }
)
task_id = response.json()["task_id"]
print(f"✓ 任务已提交: {task_id}")

# 3. 轮询任务状态
while True:
    response = requests.get(f"{BASE_URL}/tasks/{task_id}")
    task = response.json()
    status = task["status"]
    
    if status == "running":
        progress = task.get("progress", 0)
        print(f"  进度: {progress}%")
    
    if status == "completed":
        print("✓ 任务完成！")
        break
    elif status == "failed":
        print(f"✗ 任务失败: {task.get('error_message')}")
        exit(1)
    
    time.sleep(5)

# 4. 获取结果
response = requests.get(f"{BASE_URL}/results/{task_id}")
result = response.json()

print(f"\n=== 优化结果 ===")
print(f"最终能量: {result['final_energy']:.4f} eV")
print(f"能量变化: {result['energy_change']:.4f} eV")
print(f"优化步数: {result['num_steps']}")
print(f"收敛: {'是' if result['converged'] else '否'}")
print(f"体积变化: {result['volume_change_percent']:.2f}%")

# 5. 下载优化后的结构
response = requests.get(
    f"{BASE_URL}/results/{task_id}/download/structure",
    stream=True
)
with open("MOF-5-optimized.cif", "wb") as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)
print(f"\n✓ 优化结构已保存: MOF-5-optimized.cif")
```

---

### 场景 2: 多模型基准测试

```python
import requests
import concurrent.futures

BASE_URL = "http://localhost:8000/api/v1"

# 待测试的模型
models = ["mace-mp", "mace-off23", "orb-v2", "m3gnet", "chgnet"]

# 上传结构
with open("HKUST-1.cif", "r") as f:
    response = requests.post(
        f"{BASE_URL}/structures/upload",
        data={"name": "HKUST-1", "format": "cif", "content": f.read()}
    )
structure_id = response.json()["id"]

# 并行提交所有任务
def submit_task(model_key):
    response = requests.post(
        f"{BASE_URL}/tasks/optimization",
        json={
            "structure_id": structure_id,
            "model_key": model_key,
            "fmax": 0.05,
            "max_steps": 200
        }
    )
    return model_key, response.json()["task_id"]

with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
    tasks = dict(executor.map(lambda m: submit_task(m), models))

print(f"已提交 {len(tasks)} 个任务")

# 等待所有任务完成
results = {}
while tasks:
    for model_key, task_id in list(tasks.items()):
        response = requests.get(f"{BASE_URL}/tasks/{task_id}")
        task = response.json()
        
        if task["status"] == "completed":
            result_response = requests.get(f"{BASE_URL}/results/{task_id}")
            results[model_key] = result_response.json()
            del tasks[model_key]
            print(f"✓ {model_key} 完成")
        elif task["status"] == "failed":
            print(f"✗ {model_key} 失败")
            del tasks[model_key]
    
    if tasks:
        time.sleep(10)

# 输出对比结果
print("\n=== 基准测试结果 ===")
print(f"{'模型':<20} {'最终能量 (eV)':<15} {'优化步数':<10} {'用时 (s)':<10}")
print("-" * 60)
for model_key, result in results.items():
    print(f"{model_key:<20} {result['final_energy']:<15.4f} "
          f"{result['num_steps']:<10} {result['computation_time_seconds']:<10.1f}")
```

---

## 错误处理

### HTTP 状态码

| 状态码 | 说明 |
|--------|------|
| 200 | 成功 |
| 201 | 创建成功 |
| 400 | 请求参数错误 |
| 404 | 资源不存在 |
| 422 | 数据验证失败 |
| 500 | 服务器内部错误 |
| 503 | 服务不可用 |

### 错误响应格式

```json
{
  "detail": "错误详细信息",
  "error_code": "INVALID_STRUCTURE",
  "field": "structure_id"
}
```

### Python 错误处理示例

```python
import requests
from requests.exceptions import Timeout, ConnectionError

try:
    response = requests.post(
        f"{BASE_URL}/tasks/optimization",
        json=request_data,
        timeout=30
    )
    response.raise_for_status()
    task = response.json()
    
except Timeout:
    print("请求超时，请检查网络连接")
except ConnectionError:
    print("无法连接到服务器，请确认服务是否运行")
except requests.HTTPError as e:
    if e.response.status_code == 400:
        print(f"请求参数错误: {e.response.json()['detail']}")
    elif e.response.status_code == 404:
        print("结构或模型不存在")
    else:
        print(f"HTTP 错误: {e}")
except Exception as e:
    print(f"未知错误: {e}")
```

---

## 最佳实践

### 1. 超时设置

长时间计算任务需要设置合适的超时：

```python
# 优化任务：5-10 分钟
response = requests.post(url, json=data, timeout=600)

# MD 模拟：30-60 分钟
response = requests.post(url, json=data, timeout=3600)

# 声子计算：1-2 小时
response = requests.post(url, json=data, timeout=7200)
```

### 2. 任务状态轮询

建议使用指数退避策略：

```python
import time

max_wait = 60
wait_time = 5

while True:
    response = requests.get(f"{BASE_URL}/tasks/{task_id}")
    task = response.json()
    
    if task["status"] in ["completed", "failed"]:
        break
    
    time.sleep(wait_time)
    wait_time = min(wait_time * 1.5, max_wait)  # 指数增长，上限 60 秒
```

### 3. 批量任务处理

使用异步客户端提高效率：

```python
import httpx
import asyncio

async def submit_and_wait(client, structure_id, model_key):
    # 提交任务
    response = await client.post(
        f"{BASE_URL}/tasks/optimization",
        json={"structure_id": structure_id, "model_key": model_key}
    )
    task_id = response.json()["task_id"]
    
    # 等待完成
    while True:
        response = await client.get(f"{BASE_URL}/tasks/{task_id}")
        if response.json()["status"] == "completed":
            break
        await asyncio.sleep(10)
    
    # 获取结果
    response = await client.get(f"{BASE_URL}/results/{task_id}")
    return response.json()

async def main():
    async with httpx.AsyncClient(timeout=3600) as client:
        tasks = [
            submit_and_wait(client, structure_id, model)
            for model in ["mace-mp", "orb-v2", "m3gnet"]
        ]
        results = await asyncio.gather(*tasks)
        return results

results = asyncio.run(main())
```

### 4. 环境变量最佳实践

创建配置文件：

```python
# config.py
import os

class Config:
    GATEWAY_URL = os.getenv("MIRA_GATEWAY_URL", "http://localhost:8000")
    TIMEOUT = int(os.getenv("MIRA_TIMEOUT", "600"))
    RETRY_ATTEMPTS = int(os.getenv("MIRA_RETRY", "3"))
    
config = Config()

# 使用
from config import config
response = requests.post(url, json=data, timeout=config.TIMEOUT)
```

### 5. 结果缓存

避免重复计算：

```python
import json
import hashlib

def compute_hash(structure_id, model_key, params):
    data = f"{structure_id}_{model_key}_{json.dumps(params, sort_keys=True)}"
    return hashlib.md5(data.encode()).hexdigest()

# 检查缓存
cache_key = compute_hash(structure_id, model_key, optimization_params)
if cache_key in cache:
    result = cache[cache_key]
else:
    # 提交新任务
    ...
    cache[cache_key] = result
```

---

## 附录

### 支持的模型列表 (CPU 模式)

| 模型 | CPU 支持 | 说明 |
|------|----------|------|
| mace-mp | ✅ | 完全支持 |
| mace-off23 | ✅ | 完全支持 |
| orb-v2 | ✅ | 完全支持 |
| orb-v3 | ✅ | 完全支持 |
| m3gnet | ✅ | 完全支持 |
| chgnet | ✅ | 完全支持 |
| sevennet-* | ❌ | 仅 GPU |
| grace-* | ❌ | 仅 GPU |
| mattersim | ❌ | 仅 GPU |

### 环境变量完整列表

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `MIRA_GATEWAY_URL` | http://localhost:8000 | Gateway 地址 |
| `MIRA_TIMEOUT` | 600 | 请求超时（秒） |
| `MIRA_RETRY` | 3 | 重试次数 |
| `MIRA_VERBOSE` | false | 详细日志 |

---

**文档版本:** 1.0.0  
**更新日期:** 2025-12-30  
**作者:** 李世博 (Shibo Li)  
**联系方式:** shadow.li981@gmail.com
