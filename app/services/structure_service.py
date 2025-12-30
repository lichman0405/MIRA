"""
结构管理服务
"""
from typing import Optional, List, Dict, Any
from pathlib import Path
from datetime import datetime
import uuid
import json

from ase import Atoms
from ase.io import write

from ..schemas.structure import StructureInfo, StructureUpload, StructurePreview
from ..core.ase_utils import structure_from_content, get_structure_info, load_structure


class StructureService:
    """结构管理服务"""
    
    def __init__(self, structures_dir: Path):
        self.structures_dir = Path(structures_dir)
        self.structures_dir.mkdir(parents=True, exist_ok=True)
        
        # 结构元数据存储（生产环境应使用数据库）
        self._structures: Dict[str, StructureInfo] = {}
        self._metadata_file = self.structures_dir / "_metadata.json"
        
        # 加载已有结构
        self._load_metadata()
    
    def _load_metadata(self) -> None:
        """从文件加载元数据"""
        if self._metadata_file.exists():
            try:
                with open(self._metadata_file, 'r') as f:
                    data = json.load(f)
                    for key, value in data.items():
                        value['created_at'] = datetime.fromisoformat(value['created_at'])
                        self._structures[key] = StructureInfo(**value)
            except Exception:
                pass
    
    def _save_metadata(self) -> None:
        """保存元数据到文件"""
        data = {}
        for key, info in self._structures.items():
            d = info.model_dump()
            d['created_at'] = d['created_at'].isoformat()
            data[key] = d
        
        with open(self._metadata_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    async def upload_structure(self, upload: StructureUpload) -> StructureInfo:
        """
        上传结构
        
        Args:
            upload: 结构上传请求
            
        Returns:
            结构信息
        """
        # 从内容创建 Atoms 对象
        atoms = structure_from_content(
            content=upload.content,
            format=upload.format
        )
        
        # 生成 ID
        structure_id = str(uuid.uuid4())
        
        # 获取结构信息
        info = get_structure_info(atoms)
        
        # 保存结构文件
        file_path = self.structures_dir / f"{structure_id}.cif"
        write(str(file_path), atoms)
        
        # 创建结构信息
        structure_info = StructureInfo(
            id=structure_id,
            name=upload.name,
            formula=info['formula'],
            num_atoms=info['num_atoms'],
            cell_volume=info['cell_volume'],
            cell_parameters=info['cell_parameters'],
            space_group=info.get('space_group'),
            elements=info['elements'],
            file_path=str(file_path),
            metadata=upload.metadata
        )
        
        # 存储
        self._structures[structure_id] = structure_info
        self._save_metadata()
        
        return structure_info
    
    async def get_structure(self, structure_id: str) -> Optional[StructureInfo]:
        """获取结构信息"""
        return self._structures.get(structure_id)
    
    async def list_structures(
        self,
        skip: int = 0,
        limit: int = 100
    ) -> List[StructureInfo]:
        """列出所有结构"""
        structures = list(self._structures.values())
        return structures[skip:skip + limit]
    
    async def delete_structure(self, structure_id: str) -> bool:
        """删除结构"""
        if structure_id not in self._structures:
            return False
        
        info = self._structures[structure_id]
        
        # 删除文件
        file_path = Path(info.file_path)
        if file_path.exists():
            file_path.unlink()
        
        # 从存储中移除
        del self._structures[structure_id]
        self._save_metadata()
        
        return True
    
    async def get_structure_preview(
        self,
        structure_id: str
    ) -> Optional[StructurePreview]:
        """获取结构预览数据（用于可视化）"""
        info = self._structures.get(structure_id)
        if info is None:
            return None
        
        # 加载结构
        atoms = load_structure(info.file_path)
        
        return StructurePreview(
            id=structure_id,
            name=info.name,
            formula=info.formula,
            positions=atoms.get_positions().tolist(),
            symbols=atoms.get_chemical_symbols(),
            cell=atoms.get_cell().tolist(),
            pbc=atoms.pbc.tolist()
        )
    
    def get_atoms(self, structure_id: str) -> Optional[Atoms]:
        """获取 ASE Atoms 对象"""
        info = self._structures.get(structure_id)
        if info is None:
            return None
        return load_structure(info.file_path)
    
    @property
    def count(self) -> int:
        """结构数量"""
        return len(self._structures)
