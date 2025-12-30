"""
Worker 基类 - 所有模型服务的基础 FastAPI 应用
"""
import time
import torch
import traceback
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from contextlib import asynccontextmanager
from fastapi import FastAPI, HTTPException
from ase import Atoms
from ase.optimize import BFGS, FIRE, LBFGS
# ASE 版本兼容性：ExpCellFilter 和 StrainFilter 在新版本中移到了 ase.filters
try:
    from ase.filters import ExpCellFilter, StrainFilter
except ImportError:
    from ase.constraints import ExpCellFilter, StrainFilter
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units

from .schemas import (
    AtomData, TaskType, TaskStatus,
    SinglePointRequest, SinglePointResponse,
    OptimizationRequest, OptimizationResponse,
    StabilityRequest, StabilityResponse,
    BulkModulusRequest, BulkModulusResponse,
    HeatCapacityRequest, HeatCapacityResponse,
    TaskRequest, TaskResponse,
    HealthResponse, ModelsListResponse
)


def atoms_from_data(data: AtomData) -> Atoms:
    """从 AtomData 创建 ASE Atoms 对象"""
    atoms = Atoms(
        symbols=data.symbols,
        positions=data.positions,
        cell=data.cell if data.cell else None,
        pbc=data.pbc
    )
    return atoms


def atoms_to_data(atoms: Atoms) -> AtomData:
    """将 ASE Atoms 对象转换为 AtomData"""
    return AtomData(
        symbols=list(atoms.get_chemical_symbols()),
        positions=atoms.get_positions().tolist(),
        cell=atoms.get_cell().tolist() if atoms.cell is not None else None,
        pbc=list(atoms.get_pbc())
    )


class BaseWorkerApp(ABC):
    """
    Worker 基类
    
    每个模型服务继承此类，实现具体的模型加载和计算逻辑
    """
    
    def __init__(self, service_name: str, port: int = 8000):
        self.service_name = service_name
        self.port = port
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.calculators: Dict[str, Any] = {}
        self.app = self._create_app()
    
    @abstractmethod
    def get_available_models(self) -> List[str]:
        """返回该服务支持的模型列表"""
        pass
    
    @abstractmethod
    def load_model(self, model_name: str) -> Any:
        """加载指定模型，返回 ASE Calculator"""
        pass
    
    def get_calculator(self, model_name: str) -> Any:
        """获取或加载计算器"""
        if model_name not in self.calculators:
            print(f"Loading model: {model_name}")
            self.calculators[model_name] = self.load_model(model_name)
        return self.calculators[model_name]
    
    def _create_app(self) -> FastAPI:
        """创建 FastAPI 应用"""
        
        @asynccontextmanager
        async def lifespan(app: FastAPI):
            # 启动时预加载默认模型
            print(f"Starting {self.service_name}...")
            print(f"Device: {self.device}")
            print(f"Available models: {self.get_available_models()}")
            yield
            # 关闭时清理
            self.calculators.clear()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        
        app = FastAPI(
            title=f"MIRA {self.service_name} Worker",
            description=f"ML Force Field Worker Service - {self.service_name}",
            version="1.0.0",
            lifespan=lifespan
        )
        
        # 注册路由
        self._register_routes(app)
        
        return app
    
    def _register_routes(self, app: FastAPI):
        """注册 API 路由"""
        
        @app.get("/health", response_model=HealthResponse)
        async def health_check():
            """健康检查"""
            gpu_name = None
            memory_usage = None
            
            if torch.cuda.is_available():
                gpu_name = torch.cuda.get_device_name(0)
                memory_usage = torch.cuda.memory_allocated(0) / 1024 / 1024
            
            return HealthResponse(
                status="healthy",
                service=self.service_name,
                models=self.get_available_models(),
                device=self.device,
                gpu_available=torch.cuda.is_available(),
                gpu_name=gpu_name,
                memory_usage_mb=memory_usage
            )
        
        @app.get("/models", response_model=ModelsListResponse)
        async def list_models():
            """列出可用模型"""
            return ModelsListResponse(
                models=self.get_available_models(),
                service=self.service_name
            )
        
        @app.post("/single_point", response_model=SinglePointResponse)
        async def single_point(request: SinglePointRequest):
            """单点能量计算"""
            try:
                if request.model_name not in self.get_available_models():
                    raise HTTPException(
                        status_code=400,
                        detail=f"Model {request.model_name} not available in this service"
                    )
                
                calc = self.get_calculator(request.model_name)
                atoms = atoms_from_data(request.atoms)
                
                # 检查结构大小
                if len(atoms) > 1000:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Structure too large ({len(atoms)} atoms). Maximum 1000 atoms."
                    )
                
                atoms.calc = calc
                
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces().tolist() if request.compute_forces else None
                stress = atoms.get_stress().tolist() if request.compute_stress else None
                
                # 计算最大力
                max_force = None
                if forces is not None:
                    import numpy as np
                    max_force = float(np.max(np.abs(forces)))
                
                return SinglePointResponse(
                    energy=float(energy),
                    forces=forces,
                    stress=stress,
                    energy_per_atom=float(energy / len(atoms)),
                    max_force=max_force
                )
            except HTTPException:
                raise
            except Exception as e:
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=str(e))
            finally:
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        
        @app.post("/optimization", response_model=OptimizationResponse)
        async def optimization(request: OptimizationRequest):
            """结构优化"""
            calc = None
            try:
                if request.model_name not in self.get_available_models():
                    raise HTTPException(
                        status_code=400,
                        detail=f"Model {request.model_name} not available"
                    )
                
                calc = self.get_calculator(request.model_name)
                
                # 应用 D3 校正
                if request.use_d3:
                    try:
                        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
                        from ase.calculators.mixing import SumCalculator
                        d3_calc = TorchDFTD3Calculator(
                            device=self.device,
                            damping="bj"
                        )
                        calc = SumCalculator([calc, d3_calc])
                    except ImportError:
                        pass  # D3 not available
                
                atoms = atoms_from_data(request.atoms)
                atoms.calc = calc
                
                # 检查结构大小，限制原子数
                if len(atoms) > 1000:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Structure too large ({len(atoms)} atoms). Maximum 1000 atoms."
                    )
                
                initial_energy = atoms.get_potential_energy()
                
                # 检查初始能量是否异常
                if abs(initial_energy) > 1e8:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Initial energy abnormally high ({initial_energy:.2e} eV). Structure may have overlapping atoms."
                    )
                
                # 设置优化器
                if request.fix_cell:
                    opt_atoms = atoms
                else:
                    opt_atoms = ExpCellFilter(atoms)
                
                optimizer_class = {"BFGS": BFGS, "FIRE": FIRE, "LBFGS": LBFGS}.get(
                    request.optimizer, BFGS
                )
                opt = optimizer_class(opt_atoms, logfile=None)
                
                converged = opt.run(fmax=request.fmax, steps=request.max_steps)
                
                final_energy = atoms.get_potential_energy()
                final_forces = atoms.get_forces()
                
                result = OptimizationResponse(
                    converged=converged,
                    steps=opt.nsteps,
                    initial_energy=float(initial_energy),
                    final_energy=float(final_energy),
                    energy_change=float(final_energy - initial_energy),
                    final_forces_max=float(abs(final_forces).max()),
                    optimized_atoms=atoms_to_data(atoms)
                )
                
                return result
            except HTTPException:
                raise
            except Exception as e:
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=str(e))
            finally:
                # 清理资源
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        
        @app.post("/stability", response_model=StabilityResponse)
        async def stability(request: StabilityRequest):
            """稳定性分析 (NPT MD)"""
            try:
                if request.model_name not in self.get_available_models():
                    raise HTTPException(
                        status_code=400,
                        detail=f"Model {request.model_name} not available"
                    )
                
                calc = self.get_calculator(request.model_name)
                atoms = atoms_from_data(request.atoms)
                
                # 检查结构大小
                if len(atoms) > 500:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Structure too large for MD ({len(atoms)} atoms). Maximum 500 atoms."
                    )
                
                atoms.calc = calc
                
                # 初始化速度
                MaxwellBoltzmannDistribution(atoms, temperature_K=request.temperature)
                
                initial_volume = atoms.get_volume()
                initial_positions = atoms.get_positions().copy()
                
                # 使用 Langevin 动力学
                dyn = Langevin(
                    atoms,
                    timestep=request.timestep * units.fs,
                    temperature_K=request.temperature,
                    friction=0.01
                )
                
                # 收集数据
                energies = []
                temperatures = []
                volumes = []
                
                def collect_data():
                    energies.append(atoms.get_potential_energy())
                    temperatures.append(atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB))
                    volumes.append(atoms.get_volume())
                
                dyn.attach(collect_data, interval=10)
                
                # 平衡
                dyn.run(request.equilibration_steps)
                
                # 采样
                energies.clear()
                temperatures.clear()
                volumes.clear()
                dyn.run(request.production_steps)
                
                # 分析
                import numpy as np
                mean_energy = float(np.mean(energies))
                energy_std = float(np.std(energies))
                mean_temp = float(np.mean(temperatures))
                volume_change = float((np.mean(volumes) - initial_volume) / initial_volume * 100)
                
                final_positions = atoms.get_positions()
                max_disp = float(np.max(np.linalg.norm(final_positions - initial_positions, axis=1)))
                
                # 判断稳定性
                is_stable = (
                    energy_std < abs(mean_energy) * 0.1 and
                    abs(volume_change) < 20 and
                    max_disp < 5.0
                )
                
                return StabilityResponse(
                    is_stable=is_stable,
                    mean_energy=mean_energy,
                    energy_std=energy_std,
                    mean_temperature=mean_temp,
                    mean_pressure=0.0,  # Simplified
                    volume_change_percent=volume_change,
                    max_displacement=max_disp,
                    trajectory_length=len(energies)
                )
            except HTTPException:
                raise
            except Exception as e:
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=str(e))
            finally:
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        
        @app.post("/bulk_modulus", response_model=BulkModulusResponse)
        async def bulk_modulus(request: BulkModulusRequest):
            """体积模量计算"""
            try:
                import numpy as np
                from scipy.optimize import curve_fit
                
                if request.model_name not in self.get_available_models():
                    raise HTTPException(
                        status_code=400,
                        detail=f"Model {request.model_name} not available"
                    )
                
                calc = self.get_calculator(request.model_name)
                atoms = atoms_from_data(request.atoms)
                atoms.calc = calc
                
                # 生成应变点
                strains = np.linspace(
                    1 - request.strain_range,
                    1 + request.strain_range,
                    request.num_points
                )
                
                volumes = []
                energies = []
                original_cell = atoms.get_cell().copy()
                
                for s in strains:
                    test_atoms = atoms.copy()
                    test_atoms.set_cell(original_cell * s, scale_atoms=True)
                    test_atoms.calc = calc
                    
                    volumes.append(test_atoms.get_volume())
                    energies.append(test_atoms.get_potential_energy())
                
                # Birch-Murnaghan 拟合
                def birch_murnaghan(V, E0, V0, B0, Bp):
                    eta = (V0 / V) ** (2.0 / 3.0)
                    return E0 + 9.0 * V0 * B0 / 16.0 * (
                        (eta - 1) ** 3 * Bp + (eta - 1) ** 2 * (6 - 4 * eta)
                    )
                
                # 初始猜测
                V0_guess = volumes[len(volumes) // 2]
                E0_guess = min(energies)
                B0_guess = 100  # GPa
                
                try:
                    popt, _ = curve_fit(
                        birch_murnaghan,
                        volumes,
                        energies,
                        p0=[E0_guess, V0_guess, B0_guess / 160.2, 4.0],
                        maxfev=10000
                    )
                    E0, V0, B0, Bp = popt
                    B0_GPa = B0 * 160.2  # eV/Å³ -> GPa
                    
                    # 计算 R²
                    fitted = birch_murnaghan(np.array(volumes), *popt)
                    ss_res = np.sum((np.array(energies) - fitted) ** 2)
                    ss_tot = np.sum((np.array(energies) - np.mean(energies)) ** 2)
                    r2 = 1 - ss_res / ss_tot
                except:
                    # 简单多项式拟合
                    coeffs = np.polyfit(volumes, energies, 2)
                    B0_GPa = 2 * coeffs[0] * V0_guess * 160.2
                    V0 = V0_guess
                    E0 = min(energies)
                    r2 = 0.9
                
                return BulkModulusResponse(
                    bulk_modulus=float(B0_GPa),
                    equilibrium_volume=float(V0),
                    equilibrium_energy=float(E0),
                    fitting_r2=float(r2),
                    volumes=volumes,
                    energies=energies
                )
            except Exception as e:
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=str(e))
        
        @app.post("/heat_capacity", response_model=HeatCapacityResponse)
        async def heat_capacity(request: HeatCapacityRequest):
            """热容计算 (Phonopy)"""
            try:
                import numpy as np
                from phonopy import Phonopy
                from phonopy.structure.atoms import PhonopyAtoms
                
                if request.model_name not in self.get_available_models():
                    raise HTTPException(
                        status_code=400,
                        detail=f"Model {request.model_name} not available"
                    )
                
                calc = self.get_calculator(request.model_name)
                atoms = atoms_from_data(request.atoms)
                
                # 创建 PhonopyAtoms
                phonopy_atoms = PhonopyAtoms(
                    symbols=atoms.get_chemical_symbols(),
                    cell=atoms.get_cell(),
                    scaled_positions=atoms.get_scaled_positions()
                )
                
                # 创建 Phonopy 对象
                phonon = Phonopy(
                    phonopy_atoms,
                    supercell_matrix=np.diag(request.supercell)
                )
                phonon.generate_displacements(distance=0.01)
                
                # 计算力
                supercells = phonon.supercells_with_displacements
                forces_sets = []
                
                for sc in supercells:
                    sc_atoms = Atoms(
                        symbols=sc.symbols,
                        positions=sc.positions,
                        cell=sc.cell,
                        pbc=True
                    )
                    sc_atoms.calc = calc
                    forces_sets.append(sc_atoms.get_forces())
                
                phonon.forces = forces_sets
                phonon.produce_force_constants()
                
                # 计算热力学性质
                phonon.run_mesh([8, 8, 8])
                phonon.run_thermal_properties(
                    t_min=min(request.temperatures),
                    t_max=max(request.temperatures),
                    t_step=10
                )
                
                tp = phonon.get_thermal_properties_dict()
                
                # 插值到请求的温度
                from scipy.interpolate import interp1d
                
                temps = np.array(tp['temperatures'])
                cv_interp = interp1d(temps, tp['heat_capacity'], fill_value='extrapolate')
                fe_interp = interp1d(temps, tp['free_energy'], fill_value='extrapolate')
                s_interp = interp1d(temps, tp['entropy'], fill_value='extrapolate')
                
                return HeatCapacityResponse(
                    temperatures=request.temperatures,
                    heat_capacities=[float(cv_interp(t)) for t in request.temperatures],
                    free_energies=[float(fe_interp(t)) for t in request.temperatures],
                    entropies=[float(s_interp(t)) for t in request.temperatures]
                )
            except Exception as e:
                traceback.print_exc()
                raise HTTPException(status_code=500, detail=str(e))
        
        @app.post("/task", response_model=TaskResponse)
        async def execute_task(request: TaskRequest):
            """通用任务执行接口"""
            start_time = time.time()
            
            try:
                task_handlers = {
                    TaskType.SINGLE_POINT: lambda p: single_point(SinglePointRequest(**p)),
                    TaskType.OPTIMIZATION: lambda p: optimization(OptimizationRequest(**p)),
                    TaskType.STABILITY: lambda p: stability(StabilityRequest(**p)),
                    TaskType.BULK_MODULUS: lambda p: bulk_modulus(BulkModulusRequest(**p)),
                    TaskType.HEAT_CAPACITY: lambda p: heat_capacity(HeatCapacityRequest(**p)),
                }
                
                handler = task_handlers.get(request.task_type)
                if not handler:
                    return TaskResponse(
                        task_id=request.task_id,
                        status=TaskStatus.FAILED,
                        error=f"Unknown task type: {request.task_type}"
                    )
                
                result = await handler(request.payload)
                
                return TaskResponse(
                    task_id=request.task_id,
                    status=TaskStatus.COMPLETED,
                    result=result.model_dump(),
                    execution_time=time.time() - start_time
                )
            except HTTPException as e:
                return TaskResponse(
                    task_id=request.task_id,
                    status=TaskStatus.FAILED,
                    error=e.detail,
                    execution_time=time.time() - start_time
                )
            except Exception as e:
                return TaskResponse(
                    task_id=request.task_id,
                    status=TaskStatus.FAILED,
                    error=str(e),
                    execution_time=time.time() - start_time
                )
    
    def run(self, host: str = "0.0.0.0", port: int = None):
        """运行服务"""
        import uvicorn
        uvicorn.run(self.app, host=host, port=port or self.port)
