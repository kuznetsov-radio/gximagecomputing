"""I/O helpers for gximagecomputing models, EBTEL tables, and rendered map products."""

__all__ = [
    "resolve_ebtel_path",
    "load_ebtel",
    "load_ebtel_none",
    "save_h5_maps",
    "save_h5_euv_maps",
    "build_h5_from_sav",
    "estimate_hpc_center",
    "infer_center_from_execute",
    "infer_fov_from_execute",
    "load_model_hdf",
    "load_model_sav",
    "load_model_dict",
    "gx_box2id",
    "gx_voxelid",
]


def __getattr__(name):
    if name in {"resolve_ebtel_path", "load_ebtel", "load_ebtel_none"}:
        from .ebtel import resolve_ebtel_path, load_ebtel, load_ebtel_none

        return {
            "resolve_ebtel_path": resolve_ebtel_path,
            "load_ebtel": load_ebtel,
            "load_ebtel_none": load_ebtel_none,
        }[name]
    if name in {"save_h5_maps", "save_h5_euv_maps"}:
        from .maps_h5 import save_h5_maps, save_h5_euv_maps

        return {"save_h5_maps": save_h5_maps, "save_h5_euv_maps": save_h5_euv_maps}[name]
    if name in {
        "estimate_hpc_center",
        "infer_center_from_execute",
        "infer_fov_from_execute",
        "load_model_hdf",
        "load_model_sav",
        "load_model_dict",
    }:
        from .model import (
            estimate_hpc_center,
            infer_center_from_execute,
            infer_fov_from_execute,
            load_model_hdf,
            load_model_sav,
            load_model_dict,
        )

        return {
            "estimate_hpc_center": estimate_hpc_center,
            "infer_center_from_execute": infer_center_from_execute,
            "infer_fov_from_execute": infer_fov_from_execute,
            "load_model_hdf": load_model_hdf,
            "load_model_sav": load_model_sav,
            "load_model_dict": load_model_dict,
        }[name]
    if name == "build_h5_from_sav":
        from .sav_to_h5 import build_h5_from_sav as _build_h5_from_sav

        return _build_h5_from_sav
    if name in {"gx_box2id", "gx_voxelid"}:
        from .voxel_id import gx_box2id, gx_voxelid

        return {"gx_box2id": gx_box2id, "gx_voxelid": gx_voxelid}[name]
    raise AttributeError(f"module 'gximagecomputing.io' has no attribute {name!r}")
