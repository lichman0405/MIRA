# MIRA - MiQroEra Interatomic-potential Reliability Arena

A comprehensive FastAPI service for benchmarking machine learning interatomic potentials on Metal-Organic Frameworks (MOFs).

## Features

### Supported Models (20+ variants)

| Family | Models | Description |
|--------|--------|-------------|
| **MACE** | MP, OFF23, OMAT, MPA, ANI | Equivariant message passing |
| **ORB** | v2, v3, OMAT-v3-LoRA | Orbital-based descriptors |
| **OMAT24** | OMat, eqV2 variants | Meta FAIRChem models |
| **GRACE** | 2L, 2M | General-purpose MLFF |
| **MatterSim** | 5M | Materials simulation |
| **SevenNet** | 0, MF-ompa, l3i5 | Seven-body neural network |
| **PosEGNN** | IBM model | Position-enhanced GNN |
| **MatGL** | M3GNet, CHGNet | Graph neural networks |

### Computational Tasks

1. **Structure Optimization**
   - Optimizers: BFGS, FIRE, LBFGS
   - Cell filter: FrechetCellFilter for full relaxation
   - D3 dispersion correction support

2. **MD Stability Testing**
   - NVT equilibration (Langevin thermostat)
   - NPT production (NPTBerendsen)
   - Coordination number analysis
   - RMSD tracking

3. **Bulk Modulus Calculation**
   - E-V curve sampling
   - Birch-Murnaghan EOS fitting
   - Automatic R² quality assessment

4. **Heat Capacity**
   - Phonon calculations with Phonopy
   - Temperature-dependent Cv
   - Imaginary mode detection

5. **QMOF Energy Evaluation**
   - Single-point energy calculation
   - Comparison with DFT references

6. **Interaction Energy Analysis**
   - Host-guest decomposition
   - Component energy breakdown

## Installation

### Requirements

- Python >= 3.10
- CUDA >= 12.0 (for GPU acceleration)
- ASE >= 3.27.0 (includes NPT support)

### Quick Start

```bash
# Clone the repository
git clone https://github.com/your-org/MIRA.git
cd MIRA

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Run the server
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
```

### Docker Deployment

```bash
# Build and run with Docker Compose
docker-compose up -d

# View logs
docker-compose logs -f mira

# Stop
docker-compose down
```

## API Usage

### Base URL
```
http://localhost:8000/api/v1
```

### Interactive Documentation
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

### Example: Structure Optimization

```python
import requests

# Upload structure
with open("structure.cif", "r") as f:
    content = f.read()

response = requests.post(
    "http://localhost:8000/api/v1/structures/upload",
    data={
        "name": "MOF-5",
        "format": "cif",
        "content": content
    }
)
structure_id = response.json()["id"]

# Submit optimization task
response = requests.post(
    "http://localhost:8000/api/v1/tasks/optimization",
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

# Check progress
response = requests.get(f"http://localhost:8000/api/v1/tasks/{task_id}")
print(response.json())

# Get result when completed
response = requests.get(f"http://localhost:8000/api/v1/results/{task_id}")
result = response.json()
print(f"Final energy: {result['final_energy']} eV")
```

### Example: MD Stability Test

```python
response = requests.post(
    "http://localhost:8000/api/v1/tasks/stability",
    json={
        "structure_id": structure_id,
        "model_key": "mace-mp",
        "nvt_steps": 5000,
        "nvt_temperature": 300,
        "npt_steps": 10000,
        "npt_temperature": 300,
        "npt_pressure": 1.01325,
        "timestep": 1.0,
        "pre_optimize": True,
        "enable_d3": True
    }
)
```

## Configuration

Environment variables (or `.env` file):

| Variable | Default | Description |
|----------|---------|-------------|
| `DEBUG` | false | Enable debug mode |
| `DEFAULT_DEVICE` | cuda | Default compute device |
| `STRUCTURES_DIR` | ./data/structures | Structure storage |
| `RESULTS_DIR` | ./data/results | Results storage |
| `MAX_WORKERS` | 4 | Parallel task workers |
| `MAX_MD_STEPS` | 100000 | MD step limit |
| `MAX_OPT_STEPS` | 2000 | Optimization step limit |
| `CORS_ORIGINS` | ["*"] | Allowed CORS origins |

## Project Structure

```
MIRA/
├── app/
│   ├── __init__.py
│   ├── main.py              # FastAPI application
│   ├── config.py            # Configuration
│   ├── dependencies.py      # Dependency injection
│   ├── api/
│   │   └── v1/
│   │       ├── router.py    # API router
│   │       ├── models.py    # Model endpoints
│   │       ├── structures.py# Structure endpoints
│   │       ├── tasks.py     # Task endpoints
│   │       └── results.py   # Result endpoints
│   ├── schemas/
│   │   ├── model.py         # Model schemas
│   │   ├── structure.py     # Structure schemas
│   │   ├── task.py          # Task schemas
│   │   └── result.py        # Result schemas
│   ├── models/
│   │   ├── base.py          # Base adapter class
│   │   ├── registry.py      # Model registry
│   │   ├── mace_adapter.py  # MACE adapter
│   │   ├── orb_adapter.py   # ORB adapter
│   │   └── ...              # Other adapters
│   ├── services/
│   │   ├── optimization.py  # Optimization service
│   │   ├── stability.py     # Stability service
│   │   ├── bulk_modulus.py  # Bulk modulus service
│   │   ├── heat_capacity.py # Heat capacity service
│   │   ├── structure_service.py
│   │   └── task_service.py
│   └── core/
│       └── ase_utils.py     # ASE utilities
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

## Notes

- **GPU Memory**: Large models (MACE-MPA, SevenNet-l3i5) may require 16GB+ VRAM
- **ASE Version**: Requires ASE >= 3.27.0 for NPT dynamics support
- **D3 Correction**: Some tasks benefit from DFT-D3 dispersion correction
- **Phonon Calculations**: Heat capacity requires sufficient supercell size

## Author

**Shibo Li** (lishibo)  
📧 shadow.li981@gmail.com

## License

MIT License

## Citation

If you use MIRA in your research, please cite:

```bibtex
@software{mira2025,
  author={Li, Shibo},
  title={MIRA: MiQroEra Interatomic-potential Reliability Arena},
  year={2025},
  url={https://github.com/lishibo/MIRA}
}
```
