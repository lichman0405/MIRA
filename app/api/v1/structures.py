"""
结构管理 API 路由
"""
from typing import List, Optional
from fastapi import APIRouter, HTTPException, status, UploadFile, File, Form

from ...schemas.structure import (
    StructureInfo, StructureUpload, StructureListResponse, StructurePreview
)
from ...dependencies import get_structure_service


router = APIRouter(prefix="/structures", tags=["Structures"])


@router.post("/upload", response_model=StructureInfo)
async def upload_structure(
    name: str = Form(..., description="结构名称"),
    format: str = Form("cif", description="文件格式: cif/poscar/xyz"),
    content: str = Form(..., description="结构文件内容")
):
    """
    上传结构文件
    
    支持的格式：CIF, POSCAR/VASP, XYZ, extXYZ
    """
    service = get_structure_service()
    
    try:
        upload = StructureUpload(
            name=name,
            format=format,
            content=content
        )
        
        result = await service.upload_structure(upload)
        return result
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Failed to parse structure: {str(e)}"
        )


@router.post("/upload/file", response_model=StructureInfo)
async def upload_structure_file(
    file: UploadFile = File(..., description="结构文件"),
    name: Optional[str] = Form(None, description="结构名称（可选）")
):
    """通过文件上传结构"""
    service = get_structure_service()
    
    # 读取文件内容
    content = await file.read()
    content_str = content.decode('utf-8')
    
    # 确定格式
    filename = file.filename or "structure.cif"
    ext = filename.split('.')[-1].lower()
    format_map = {
        'cif': 'cif',
        'vasp': 'poscar',
        'poscar': 'poscar',
        'xyz': 'xyz',
        'extxyz': 'extxyz',
        'json': 'json'
    }
    format = format_map.get(ext, 'cif')
    
    try:
        upload = StructureUpload(
            name=name or filename,
            format=format,
            content=content_str
        )
        
        result = await service.upload_structure(upload)
        return result
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Failed to parse structure: {str(e)}"
        )


@router.get("/", response_model=StructureListResponse)
async def list_structures(
    skip: int = 0,
    limit: int = 100
):
    """列出所有已上传的结构"""
    service = get_structure_service()
    
    structures = await service.list_structures(skip=skip, limit=limit)
    
    return StructureListResponse(
        total=service.count,
        structures=structures
    )


@router.get("/{structure_id}", response_model=StructureInfo)
async def get_structure(structure_id: str):
    """获取结构详细信息"""
    service = get_structure_service()
    
    structure = await service.get_structure(structure_id)
    if structure is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Structure '{structure_id}' not found"
        )
    
    return structure


@router.get("/{structure_id}/preview", response_model=StructurePreview)
async def get_structure_preview(structure_id: str):
    """
    获取结构预览数据
    
    返回坐标、元素符号和晶胞信息，用于前端可视化
    """
    service = get_structure_service()
    
    preview = await service.get_structure_preview(structure_id)
    if preview is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Structure '{structure_id}' not found"
        )
    
    return preview


@router.delete("/{structure_id}")
async def delete_structure(structure_id: str):
    """删除结构"""
    service = get_structure_service()
    
    success = await service.delete_structure(structure_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Structure '{structure_id}' not found"
        )
    
    return {"status": "deleted", "structure_id": structure_id}
