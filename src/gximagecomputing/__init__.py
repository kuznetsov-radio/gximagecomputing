__all__ = [
    "GXRadioImageComputing",
    "resolve_ebtel_path",
    "estimate_hpc_center",
    "infer_fov_from_execute",
    "save_h5_maps",
    "build_h5_from_sav",
]


def __getattr__(name):
    if name == "GXRadioImageComputing":
        from .radio import GXRadioImageComputing as _GXRadioImageComputing

        return _GXRadioImageComputing

    if name in {
        "resolve_ebtel_path",
        "estimate_hpc_center",
        "infer_fov_from_execute",
        "save_h5_maps",
        "build_h5_from_sav",
    }:
        from . import io as _io

        return getattr(_io, name)

    raise AttributeError(f"module 'gximagecomputing' has no attribute {name!r}")
