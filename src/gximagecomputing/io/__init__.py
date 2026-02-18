"""I/O helpers for gximagecomputing models, EBTEL tables, and rendered map products."""

__all__ = [
    "resolve_ebtel_path",
    "save_h5_maps",
    "build_h5_from_sav",
    "estimate_hpc_center",
    "infer_fov_from_execute",
    "gx_box2id",
    "gx_voxelid",
]


def __getattr__(name):
    if name == "resolve_ebtel_path":
        from .ebtel import resolve_ebtel_path as _resolve_ebtel_path

        return _resolve_ebtel_path
    if name == "save_h5_maps":
        from .maps_h5 import save_h5_maps as _save_h5_maps

        return _save_h5_maps
    if name in {"estimate_hpc_center", "infer_fov_from_execute"}:
        from .model import estimate_hpc_center, infer_fov_from_execute

        return {"estimate_hpc_center": estimate_hpc_center, "infer_fov_from_execute": infer_fov_from_execute}[name]
    if name == "build_h5_from_sav":
        from .sav_to_h5 import build_h5_from_sav as _build_h5_from_sav

        return _build_h5_from_sav
    if name in {"gx_box2id", "gx_voxelid"}:
        from .voxel_id import gx_box2id, gx_voxelid

        return {"gx_box2id": gx_box2id, "gx_voxelid": gx_voxelid}[name]
    raise AttributeError(f"module 'gximagecomputing.io' has no attribute {name!r}")
