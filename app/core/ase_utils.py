"""
ASE 工具函数
提供结构读写、分析等通用功能
"""
from typing import Optional, List, Dict, Any, Union, Tuple
from pathlib import Path
import numpy as np
import uuid
import json

try:
    from ase import Atoms
    from ase.io import read, write
    from ase.geometry import get_duplicate_atoms
    from ase.spacegroup import get_spacegroup
except ImportError:
    raise ImportError("ASE is required. Install with: pip install ase>=3.27.0")


def load_structure(
    file_path: Union[str, Path],
    format: Optional[str] = None
) -> Atoms:
    """
    从文件加载结构
    
    Args:
        file_path: 文件路径
        format: 文件格式（可选，自动检测）
        
    Returns:
        ASE Atoms 对象
    """
    return read(str(file_path), format=format)


def save_structure(
    atoms: Atoms,
    file_path: Union[str, Path],
    format: Optional[str] = None
) -> None:
    """
    保存结构到文件
    
    Args:
        atoms: ASE Atoms 对象
        file_path: 保存路径
        format: 文件格式（可选，根据扩展名自动检测）
    """
    write(str(file_path), atoms, format=format)


def structure_from_content(
    content: str,
    format: str,
    filename: Optional[str] = None
) -> Atoms:
    """
    从字符串内容创建结构
    
    Args:
        content: 文件内容字符串
        format: 文件格式 (cif/poscar/xyz/json)
        filename: 临时文件名（可选）
        
    Returns:
        ASE Atoms 对象
    """
    import tempfile
    import os
    
    # 确定文件扩展名
    ext_map = {
        "cif": ".cif",
        "poscar": ".vasp",
        "xyz": ".xyz",
        "json": ".json",
        "extxyz": ".extxyz"
    }
    ext = ext_map.get(format.lower(), f".{format}")
    
    # 写入临时文件
    with tempfile.NamedTemporaryFile(
        mode='w', 
        suffix=ext, 
        delete=False
    ) as f:
        f.write(content)
        temp_path = f.name
    
    try:
        atoms = read(temp_path, format=format if format != "poscar" else "vasp")
    finally:
        os.unlink(temp_path)
    
    return atoms


def get_structure_info(atoms: Atoms) -> Dict[str, Any]:
    """
    获取结构的详细信息
    
    Args:
        atoms: ASE Atoms 对象
        
    Returns:
        包含结构信息的字典
    """
    # 获取晶胞参数
    cell = atoms.get_cell()
    cell_lengths = cell.lengths()
    cell_angles = cell.angles()
    
    # 获取化学式
    formula = atoms.get_chemical_formula(mode='hill')
    
    # 获取元素列表
    elements = list(set(atoms.get_chemical_symbols()))
    
    # 尝试获取空间群
    try:
        spg = get_spacegroup(atoms)
        space_group = f"{spg.symbol} ({spg.no})"
    except Exception:
        space_group = None
    
    return {
        "formula": formula,
        "num_atoms": len(atoms),
        "cell_volume": atoms.get_volume(),
        "cell_parameters": {
            "a": float(cell_lengths[0]),
            "b": float(cell_lengths[1]),
            "c": float(cell_lengths[2]),
            "alpha": float(cell_angles[0]),
            "beta": float(cell_angles[1]),
            "gamma": float(cell_angles[2])
        },
        "space_group": space_group,
        "elements": sorted(elements),
        "pbc": atoms.pbc.tolist()
    }


def calculate_rmsd(atoms1: Atoms, atoms2: Atoms) -> float:
    """
    计算两个结构之间的 RMSD
    
    Args:
        atoms1: 参考结构
        atoms2: 比较结构
        
    Returns:
        RMSD 值 (Å)
    """
    if len(atoms1) != len(atoms2):
        raise ValueError("Structures have different number of atoms")
    
    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()
    
    diff = pos1 - pos2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return float(rmsd)


def get_coordination_numbers(
    atoms: Atoms,
    metal_symbols: Optional[List[str]] = None,
    cutoff: float = 3.0
) -> Dict[str, List[int]]:
    """
    计算金属原子的配位数
    
    Args:
        atoms: ASE Atoms 对象
        metal_symbols: 金属元素符号列表（默认自动检测）
        cutoff: 配位半径截断 (Å)
        
    Returns:
        每种金属的配位数列表
    """
    from ase.neighborlist import neighbor_list
    
    # 默认金属元素
    if metal_symbols is None:
        metals = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
                  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                  'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                  'La', 'Ce', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                  'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb', 'Bi'}
        symbols = atoms.get_chemical_symbols()
        metal_symbols = [s for s in set(symbols) if s in metals]
    
    if not metal_symbols:
        return {}
    
    # 计算邻居列表
    i_list, j_list = neighbor_list('ij', atoms, cutoff)
    
    # 统计配位数
    coordination = {m: [] for m in metal_symbols}
    symbols = atoms.get_chemical_symbols()
    
    for metal in metal_symbols:
        metal_indices = [i for i, s in enumerate(symbols) if s == metal]
        for idx in metal_indices:
            cn = np.sum(i_list == idx)
            coordination[metal].append(int(cn))
    
    return coordination


def scale_cell(atoms: Atoms, scale_factor: float) -> Atoms:
    """
    等比例缩放晶胞
    
    Args:
        atoms: ASE Atoms 对象
        scale_factor: 缩放因子（体积的 1/3 次方）
        
    Returns:
        缩放后的新 Atoms 对象
    """
    scaled = atoms.copy()
    scaled.set_cell(atoms.get_cell() * scale_factor, scale_atoms=True)
    return scaled


def get_volume_scaled_structures(
    atoms: Atoms,
    volume_range: float = 0.1,
    num_points: int = 11
) -> List[Tuple[float, Atoms]]:
    """
    生成一系列不同体积的结构（用于 E-V 曲线）
    
    Args:
        atoms: 原始结构
        volume_range: 体积变化范围 (±比例)
        num_points: 采样点数
        
    Returns:
        (体积, Atoms) 元组列表
    """
    v0 = atoms.get_volume()
    scale_factors = np.linspace(1 - volume_range, 1 + volume_range, num_points)
    
    structures = []
    for sf in scale_factors:
        # 体积缩放因子的立方根
        cell_sf = sf ** (1/3)
        scaled = scale_cell(atoms, cell_sf)
        structures.append((scaled.get_volume(), scaled))
    
    return structures


def create_supercell(atoms: Atoms, supercell: List[int]) -> Atoms:
    """
    创建超胞
    
    Args:
        atoms: 原始结构
        supercell: 超胞大小 [nx, ny, nz]
        
    Returns:
        超胞结构
    """
    from ase.build import make_supercell
    import numpy as np
    
    P = np.diag(supercell)
    return make_supercell(atoms, P)
