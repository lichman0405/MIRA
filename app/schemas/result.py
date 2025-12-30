"""
结果相关的 Pydantic 数据模型
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class OptimizationResult(BaseModel):
    """结构优化结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 能量信息
    initial_energy: float = Field(..., description="初始能量 (eV)")
    final_energy: float = Field(..., description="最终能量 (eV)")
    energy_change: float = Field(..., description="能量变化 (eV)")
    energy_per_atom: float = Field(..., description="每原子能量 (eV/atom)")
    
    # 结构变化
    rmsd: float = Field(..., description="RMSD (Å)")
    initial_volume: float = Field(..., description="初始体积 (Å³)")
    final_volume: float = Field(..., description="最终体积 (Å³)")
    volume_change_percent: float = Field(..., description="体积变化百分比 (%)")
    
    # 优化信息
    num_steps: int = Field(..., description="优化步数")
    converged: bool = Field(..., description="是否收敛")
    final_fmax: float = Field(..., description="最终最大力 (eV/Å)")
    
    # 晶胞参数
    cell_parameters: Dict[str, float] = Field(..., description="晶胞参数")
    
    # 文件路径
    optimized_structure_file: str = Field(..., description="优化后结构文件路径")
    trajectory_file: Optional[str] = Field(None, description="轨迹文件路径")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "task_id": "550e8400-e29b-41d4-a716-446655440000",
                "initial_energy": -1500.234,
                "final_energy": -1502.567,
                "energy_change": -2.333,
                "energy_per_atom": -3.546,
                "rmsd": 0.234,
                "initial_volume": 17237.5,
                "final_volume": 17150.3,
                "volume_change_percent": -0.51,
                "num_steps": 45,
                "converged": True,
                "final_fmax": 0.0089,
                "cell_parameters": {
                    "a": 25.75, "b": 25.75, "c": 25.75,
                    "alpha": 90.0, "beta": 90.0, "gamma": 90.0
                },
                "optimized_structure_file": "./results/xxx_optimized.cif"
            }
        }
    }


class StabilityResult(BaseModel):
    """MD 稳定性模拟结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 温度统计
    average_temperature: float = Field(..., description="平均温度 (K)")
    temperature_std: float = Field(..., description="温度标准差 (K)")
    
    # 压力统计
    average_pressure: float = Field(..., description="平均压力 (bar)")
    pressure_std: float = Field(..., description="压力标准差 (bar)")
    
    # 结构稳定性
    final_rmsd: float = Field(..., description="最终 RMSD (Å)")
    max_rmsd: float = Field(..., description="最大 RMSD (Å)")
    volume_drift_percent: float = Field(..., description="体积漂移百分比 (%)")
    
    # 配位分析
    coordination_analysis: Dict[str, Any] = Field(..., description="配位数分析")
    
    # 能量统计
    average_energy: float = Field(..., description="平均能量 (eV)")
    energy_std: float = Field(..., description="能量标准差 (eV)")
    
    # 模拟信息
    total_steps: int = Field(..., description="总步数")
    simulation_time_ps: float = Field(..., description="模拟时间 (ps)")
    
    # 文件路径
    trajectory_file: str = Field(..., description="轨迹文件路径")
    final_structure_file: str = Field(..., description="最终结构文件路径")
    energy_plot: Optional[str] = Field(None, description="能量图路径")
    rmsd_plot: Optional[str] = Field(None, description="RMSD 图路径")
    temperature_plot: Optional[str] = Field(None, description="温度图路径")


class BulkModulusResult(BaseModel):
    """体积模量计算结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 主要结果
    bulk_modulus: float = Field(..., description="体积模量 B₀ (GPa)")
    bulk_modulus_derivative: float = Field(..., description="体积模量导数 B₀'")
    equilibrium_volume: float = Field(..., description="平衡体积 V₀ (Å³)")
    equilibrium_energy: float = Field(..., description="平衡能量 E₀ (eV)")
    
    # 拟合质量
    r_squared: float = Field(..., description="拟合 R² 值")
    residual: float = Field(..., description="拟合残差")
    
    # 数据点
    ev_curve_data: List[Dict[str, float]] = Field(..., description="E-V 曲线数据点")
    
    # 文件路径
    ev_curve_plot: Optional[str] = Field(None, description="E-V 曲线图路径")
    data_file: Optional[str] = Field(None, description="数据文件路径")
    
    model_config = {
        "json_schema_extra": {
            "example": {
                "task_id": "550e8400-e29b-41d4-a716-446655440000",
                "bulk_modulus": 15.6,
                "bulk_modulus_derivative": 4.2,
                "equilibrium_volume": 17150.3,
                "equilibrium_energy": -1502.567,
                "r_squared": 0.9998,
                "residual": 0.0012,
                "ev_curve_data": [
                    {"volume": 15435.27, "energy": -1498.234},
                    {"volume": 17150.30, "energy": -1502.567},
                    {"volume": 18865.33, "energy": -1499.123}
                ]
            }
        }
    }


class HeatCapacityResult(BaseModel):
    """热容计算结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 300K 下的热力学性质
    cv_300k: float = Field(..., description="300K 等容热容 Cv (J/(mol·K))")
    entropy_300k: float = Field(..., description="300K 熵 S (J/(mol·K))")
    free_energy_300k: float = Field(..., description="300K 自由能 F (kJ/mol)")
    internal_energy_300k: float = Field(..., description="300K 内能 U (kJ/mol)")
    
    # 温度依赖数据
    temperature_cv_data: List[Dict[str, float]] = Field(..., description="T-Cv 数据")
    
    # 声子信息
    phonon_frequencies: List[float] = Field(..., description="声子频率 (THz)")
    has_imaginary_modes: bool = Field(..., description="是否存在虚频")
    num_imaginary_modes: int = Field(0, description="虚频数量")
    
    # 文件路径
    phonon_dos_plot: Optional[str] = Field(None, description="声子 DOS 图路径")
    cv_plot: Optional[str] = Field(None, description="Cv-T 曲线图路径")
    band_structure_plot: Optional[str] = Field(None, description="声子能带图路径")
    data_file: Optional[str] = Field(None, description="数据文件路径")


class QMOFEnergyResult(BaseModel):
    """QMOF 能量对比结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 统计指标
    mae: float = Field(..., description="平均绝对误差 MAE (eV/atom)")
    rmse: float = Field(..., description="均方根误差 RMSE (eV/atom)")
    max_error: float = Field(..., description="最大误差 (eV/atom)")
    r_squared: float = Field(..., description="相关系数 R²")
    
    # 详细数据
    energy_comparison: List[Dict[str, Any]] = Field(..., description="能量对比数据")
    
    # 文件路径
    parity_plot: Optional[str] = Field(None, description="Parity 图路径")
    data_file: Optional[str] = Field(None, description="数据文件路径")


class InteractionEnergyResult(BaseModel):
    """相互作用能计算结果"""
    task_id: str = Field(..., description="任务 ID")
    
    # 统计结果
    average_interaction_energy: float = Field(..., description="平均相互作用能 (eV)")
    min_interaction_energy: float = Field(..., description="最小相互作用能 (eV)")
    max_interaction_energy: float = Field(..., description="最大相互作用能 (eV)")
    
    # 详细数据
    binding_sites_data: List[Dict[str, Any]] = Field(..., description="各位点数据")
    
    # 参考值比较（如果有）
    reference_energy: Optional[float] = Field(None, description="参考值 (eV)")
    mae: Optional[float] = Field(None, description="MAE (eV)")
    
    # 文件路径
    visualization_file: Optional[str] = Field(None, description="可视化文件路径")
    data_file: Optional[str] = Field(None, description="数据文件路径")
