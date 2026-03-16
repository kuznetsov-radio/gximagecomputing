from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np

from .workflows import render_euv as _render_euv_workflow
from .workflows import render_mw as _render_mw_workflow


@dataclass(slots=True)
class ObserverOverrides:
    """Optional observer metadata overrides applied before rendering."""

    dsun_cm: float | None = None
    lonc_deg: float | None = None
    b0sun_deg: float | None = None


@dataclass(slots=True)
class MapGeometry:
    """Optional output map geometry controls."""

    xc: float | None = None
    yc: float | None = None
    dx: float | None = None
    dy: float | None = None
    pixel_scale_arcsec: float | None = None
    nx: int | None = None
    ny: int | None = None
    xrange: tuple[float, float] | None = None
    yrange: tuple[float, float] | None = None


@dataclass(slots=True)
class CoronalPlasmaParameters:
    """Coronal plasma/heating parameters forwarded to the render engine."""

    tbase: float | None = None
    nbase: float | None = None
    q0: float | None = None
    a: float | None = None
    b: float | None = None
    mode: int = 0
    selective_heating: bool = False
    shtable: Sequence[Sequence[float]] | np.ndarray | None = None


@dataclass(slots=True)
class MWRenderOptions:
    model_path: str | Path
    model_format: Literal["h5", "sav", "auto"] = "auto"
    ebtel_path: str | None = None
    output_dir: str | Path | None = None
    output_name: str | None = None
    output_format: Literal["h5", "fits", "both"] = "h5"
    freqlist_ghz: Sequence[float] | None = None
    plasma: CoronalPlasmaParameters | None = None
    omp_threads: int = 8
    geometry: MapGeometry | None = None
    observer: ObserverOverrides | None = None
    save_outputs: bool = True
    write_preview: bool = True
    verbose: bool = False


@dataclass(slots=True)
class EUVRenderOptions:
    model_path: str | Path
    model_format: Literal["h5", "sav", "auto"] = "auto"
    ebtel_path: str | None = None
    output_dir: str | Path | None = None
    output_name: str | None = None
    channels: Sequence[str] | None = None
    instrument: str | None = None
    response_sav: str | Path | None = None
    response: Any | None = None
    response_dt: Any | None = None
    response_meta: Any | None = None
    plasma: CoronalPlasmaParameters | None = None
    omp_threads: int = 8
    geometry: MapGeometry | None = None
    observer: ObserverOverrides | None = None
    save_outputs: bool = True
    write_preview: bool = True
    verbose: bool = False


@dataclass(slots=True)
class RenderGeometryInfo:
    xc_arcsec: float
    yc_arcsec: float
    dx_arcsec: float
    dy_arcsec: float
    nx: int
    ny: int
    fov_x_arcsec: float
    fov_y_arcsec: float


@dataclass(slots=True)
class MWOutputFiles:
    output_dir: Path
    save_outputs: bool
    write_preview: bool
    h5_path: Path | None
    preview_png: Path | None
    fits_paths: list[Path]


@dataclass(slots=True)
class EUVOutputFiles:
    output_dir: Path
    save_outputs: bool
    write_preview: bool
    h5_path: Path | None
    preview_png: Path | None


@dataclass(slots=True)
class EUVResponseInfo:
    instrument: str
    channels: list[str]
    source: str
    mode: str


@dataclass(slots=True)
class MWRenderResult:
    library_path: str
    model_path: Path
    model_format: str
    ebtel_path: str
    observer_overrides_applied: dict[str, float]
    center_source: str
    geometry: RenderGeometryInfo
    obs_time_iso: str
    freqlist_ghz: list[float]
    plasma: CoronalPlasmaParameters
    ti: np.ndarray
    tv: np.ndarray
    outputs: MWOutputFiles
    raw_result: dict[str, Any]


@dataclass(slots=True)
class EUVRenderResult:
    library_path: str
    model_path: Path
    model_format: str
    ebtel_path: str
    observer_overrides_applied: dict[str, float]
    center_source: str
    geometry: RenderGeometryInfo
    obs_time_iso: str
    response: EUVResponseInfo
    plasma: CoronalPlasmaParameters
    flux_corona: np.ndarray
    flux_tr: np.ndarray
    outputs: EUVOutputFiles
    raw_result: dict[str, Any]


def _geometry_to_kwargs(geometry: MapGeometry | None) -> dict:
    g = geometry or MapGeometry()
    return {
        "xc": g.xc,
        "yc": g.yc,
        "dx": g.dx,
        "dy": g.dy,
        "pixel_scale_arcsec": g.pixel_scale_arcsec,
        "nx": g.nx,
        "ny": g.ny,
        "xrange": list(g.xrange) if g.xrange is not None else None,
        "yrange": list(g.yrange) if g.yrange is not None else None,
    }


def _observer_to_kwargs(observer: ObserverOverrides | None) -> dict:
    o = observer or ObserverOverrides()
    return {
        "dsun_cm": o.dsun_cm,
        "lonc_deg": o.lonc_deg,
        "b0sun_deg": o.b0sun_deg,
    }


def _plasma_to_kwargs(plasma: CoronalPlasmaParameters | None) -> dict:
    if plasma is None:
        return {
            "tbase": None,
            "nbase": None,
            "q0": None,
            "a": None,
            "b": None,
            "corona_mode": 0,
            "selective_heating": False,
            "shtable": None,
            "shtable_path": None,
            "force_isothermal": False,
            "interpol_b": False,
            "analytical_nt": False,
        }

    p = plasma
    return {
        "tbase": p.tbase,
        "nbase": p.nbase,
        "q0": p.q0,
        "a": p.a,
        "b": p.b,
        "corona_mode": p.mode,
        "selective_heating": bool(p.selective_heating),
        "shtable": (np.asarray(p.shtable, dtype=np.float64) if p.shtable is not None else None),
        "shtable_path": None,
        "force_isothermal": False,
        "interpol_b": False,
        "analytical_nt": False,
    }


def _plasma_from_dict(d: dict[str, Any]) -> CoronalPlasmaParameters:
    return CoronalPlasmaParameters(
        tbase=float(d["tbase_k"]),
        nbase=float(d["nbase_cm3"]),
        q0=float(d["q0"]),
        a=float(d["a"]),
        b=float(d["b"]),
        mode=int(d["mode"]),
        selective_heating=bool(d.get("selective_heating", d.get("shtable") is not None)),
        shtable=(
            None
            if d.get("shtable") is None
            else np.asarray(d["shtable"], dtype=np.float64)
        ),
    )


def _geometry_from_dict(d: dict[str, Any]) -> RenderGeometryInfo:
    return RenderGeometryInfo(
        xc_arcsec=float(d["xc_arcsec"]),
        yc_arcsec=float(d["yc_arcsec"]),
        dx_arcsec=float(d["dx_arcsec"]),
        dy_arcsec=float(d["dy_arcsec"]),
        nx=int(d["nx"]),
        ny=int(d["ny"]),
        fov_x_arcsec=float(d["fov_x_arcsec"]),
        fov_y_arcsec=float(d["fov_y_arcsec"]),
    )


def _mw_result_from_workflow(d: dict[str, Any]) -> MWRenderResult:
    outputs = d["outputs"]
    raw = d["result"]
    return MWRenderResult(
        library_path=str(d["library_path"]),
        model_path=Path(d["model_path"]),
        model_format=str(d["model_format"]),
        ebtel_path=str(d["ebtel_path"]),
        observer_overrides_applied=dict(d.get("observer_overrides_applied", {})),
        center_source=str(d["center_source"]),
        geometry=_geometry_from_dict(d["geometry"]),
        obs_time_iso=str(d["obs_time_iso"]),
        freqlist_ghz=[float(x) for x in d["freqlist_ghz"]],
        plasma=_plasma_from_dict(d["plasma"]),
        ti=np.asarray(raw["TI"]),
        tv=np.asarray(raw["TV"]),
        outputs=MWOutputFiles(
            output_dir=Path(outputs["output_dir"]),
            save_outputs=bool(outputs.get("save_outputs", True)),
            write_preview=bool(outputs.get("write_preview", True)),
            h5_path=(Path(outputs["h5_path"]) if outputs.get("h5_path") else None),
            preview_png=(Path(outputs["preview_png"]) if outputs.get("preview_png") else None),
            fits_paths=[Path(p) for p in outputs.get("fits_paths", [])],
        ),
        raw_result=raw,
    )


def _euv_result_from_workflow(d: dict[str, Any]) -> EUVRenderResult:
    outputs = d["outputs"]
    raw = d["result"]
    resp = d["response"]
    return EUVRenderResult(
        library_path=str(d["library_path"]),
        model_path=Path(d["model_path"]),
        model_format=str(d["model_format"]),
        ebtel_path=str(d["ebtel_path"]),
        observer_overrides_applied=dict(d.get("observer_overrides_applied", {})),
        center_source=str(d["center_source"]),
        geometry=_geometry_from_dict(d["geometry"]),
        obs_time_iso=str(d["obs_time_iso"]),
        response=EUVResponseInfo(
            instrument=str(resp["instrument"]),
            channels=[str(c) for c in resp["channels"]],
            source=str(resp["source"]),
            mode=str(resp["mode"]),
        ),
        plasma=_plasma_from_dict(d["plasma"]),
        flux_corona=np.asarray(raw["flux_corona"]),
        flux_tr=np.asarray(raw["flux_tr"]),
        outputs=EUVOutputFiles(
            output_dir=Path(outputs["output_dir"]),
            save_outputs=bool(outputs.get("save_outputs", True)),
            write_preview=bool(outputs.get("write_preview", True)),
            h5_path=(Path(outputs["h5_path"]) if outputs.get("h5_path") else None),
            preview_png=(Path(outputs["preview_png"]) if outputs.get("preview_png") else None),
        ),
        raw_result=raw,
    )


def render_mw_maps(options: MWRenderOptions) -> MWRenderResult:
    """Render microwave maps programmatically without argparse/CLI overhead.

    Returns a typed `MWRenderResult`.
    """

    ns = Namespace(
        model_path=Path(options.model_path),
        model_format=str(options.model_format),
        ebtel_path=options.ebtel_path,
        output_dir=(Path(options.output_dir) if options.output_dir is not None else _render_mw_workflow.DEFAULT_OUTDIR),
        output_name=options.output_name,
        output_format=str(options.output_format),
        frequencies_ghz=(list(options.freqlist_ghz) if options.freqlist_ghz is not None else None),
        omp_threads=int(options.omp_threads),
        save_outputs=bool(options.save_outputs),
        write_preview=bool(options.write_preview),
        **_geometry_to_kwargs(options.geometry),
        **_observer_to_kwargs(options.observer),
        **_plasma_to_kwargs(options.plasma),
    )
    return _mw_result_from_workflow(_render_mw_workflow.run(ns, verbose=bool(options.verbose)))


def render_euv_maps(options: EUVRenderOptions) -> EUVRenderResult:
    """Render EUV maps programmatically without argparse/CLI overhead.

    Returns a typed `EUVRenderResult`.
    """

    ns = Namespace(
        model_path=Path(options.model_path),
        model_format=str(options.model_format),
        ebtel_path=options.ebtel_path,
        output_dir=(Path(options.output_dir) if options.output_dir is not None else _render_euv_workflow.DEFAULT_OUTDIR),
        output_name=options.output_name,
        channels=(list(options.channels) if options.channels is not None else None),
        instrument=(str(options.instrument) if options.instrument is not None else None),
        response_sav=(Path(options.response_sav) if options.response_sav is not None else None),
        response=options.response,
        response_dt=options.response_dt,
        response_meta=options.response_meta,
        omp_threads=int(options.omp_threads),
        save_outputs=bool(options.save_outputs),
        write_preview=bool(options.write_preview),
        **_geometry_to_kwargs(options.geometry),
        **_observer_to_kwargs(options.observer),
        **_plasma_to_kwargs(options.plasma),
    )
    return _euv_result_from_workflow(_render_euv_workflow.run(ns, verbose=bool(options.verbose)))


__all__ = [
    "ObserverOverrides",
    "MapGeometry",
    "CoronalPlasmaParameters",
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
