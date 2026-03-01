__all__ = [
    "GXRadioImageComputing",
    "GXEUVImageComputing",
    "resolve_ebtel_path",
    "estimate_hpc_center",
    "infer_fov_from_execute",
    "save_h5_maps",
    "save_h5_euv_maps",
    "build_h5_from_sav",
    "ObserverOverrides",
    "MapGeometry",
    "MWRenderOptions",
    "EUVRenderOptions",
    "RenderGeometryInfo",
    "MWOutputFiles",
    "EUVOutputFiles",
    "EUVResponseInfo",
    "MWRenderResult",
    "EUVRenderResult",
    "render_mw_maps",
    "render_euv_maps",
]


def __getattr__(name):
    if name == "GXRadioImageComputing":
        from .radio import GXRadioImageComputing as _GXRadioImageComputing

        return _GXRadioImageComputing
    if name == "GXEUVImageComputing":
        from .euv import GXEUVImageComputing as _GXEUVImageComputing

        return _GXEUVImageComputing

    if name in {
        "resolve_ebtel_path",
        "estimate_hpc_center",
        "infer_fov_from_execute",
        "save_h5_maps",
        "save_h5_euv_maps",
        "build_h5_from_sav",
    }:
        from . import io as _io

        return getattr(_io, name)

    if name in {
        "ObserverOverrides",
        "MapGeometry",
        "MWRenderOptions",
        "EUVRenderOptions",
        "RenderGeometryInfo",
        "MWOutputFiles",
        "EUVOutputFiles",
        "EUVResponseInfo",
        "MWRenderResult",
        "EUVRenderResult",
        "render_mw_maps",
        "render_euv_maps",
    }:
        from . import sdk as _sdk

        return getattr(_sdk, name)

    raise AttributeError(f"module 'gximagecomputing' has no attribute {name!r}")
