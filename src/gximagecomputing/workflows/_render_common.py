from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.time import Time

from gximagecomputing.io.ebtel import load_ebtel, load_ebtel_none, resolve_ebtel_path
from gximagecomputing.io.model import (
    estimate_hpc_center,
    infer_center_from_execute,
    infer_fov_from_execute,
    load_model_hdf,
    load_model_sav,
)

DEFAULT_OUTDIR = (
    Path("C:/Temp/gximagecomputing_validation_groundtruth")
    if os.name == "nt"
    else Path("/tmp/gximagecomputing_validation_groundtruth")
)


@dataclass
class CommonRenderInputs:
    model_path: Path
    loader: str
    model: Any
    model_dt: Any
    ebtel_path: str
    ebtel_c: Any
    ebtel_dt: Any
    center_source: str
    xc: float
    yc: float
    dx: float
    dy: float
    nx: int
    ny: int
    fov_x: float
    fov_y: float
    observer_overrides_applied: dict[str, float]


@dataclass
class PlasmaDefaults:
    tbase: float
    nbase: float
    q0: float
    a: float
    b: float
    shtable: np.ndarray


def resolve_model_loader(model_path: Path, model_format: str) -> str:
    if model_format == "auto":
        return "h5" if model_path.suffix.lower() in {".h5", ".hdf5"} else "sav"
    return str(model_format)


def load_ebtel_from_arg(ebtel_path_arg: str | None) -> tuple[str, Any, Any]:
    ebtel_arg = None
    disable_ebtel = False
    if ebtel_path_arg is not None:
        ebtel_raw = str(ebtel_path_arg).strip()
        if ebtel_raw in {"", '""', "''"}:
            disable_ebtel = True
        else:
            ebtel_arg = Path(ebtel_raw)

    ebtel_env = os.environ.get("GXIMAGECOMPUTING_EBTEL_PATH", "").strip()
    if disable_ebtel or (ebtel_arg is None and not ebtel_env):
        ebtel_path = ""
        ebtel_c, ebtel_dt = load_ebtel_none()
    else:
        ebtel_path = str(resolve_ebtel_path(ebtel_arg))
        ebtel_c, ebtel_dt = load_ebtel(ebtel_path)
    return ebtel_path, ebtel_c, ebtel_dt


def load_model_and_fov(
    model_path: Path,
    loader: str,
    *,
    prefer_execute_center: bool = True,
    observer_overrides: dict[str, float | None] | None = None,
) -> tuple[Any, Any, str, float, float, float, float, dict[str, float]]:
    if loader == "h5":
        model, model_dt = load_model_hdf(str(model_path))
    else:
        model, model_dt = load_model_sav(str(model_path))
    applied_overrides: dict[str, float] = {}
    if observer_overrides:
        applied_overrides = apply_model_observer_overrides(
            model,
            dsun_cm=observer_overrides.get("dsun_cm"),
            lonc_deg=observer_overrides.get("lonc_deg"),
            b0sun_deg=observer_overrides.get("b0sun_deg"),
        )

    center_exec = infer_center_from_execute(loader_name=loader, model_path=model_path) if prefer_execute_center else None
    if center_exec is not None:
        xc_auto, yc_auto = float(center_exec[0]), float(center_exec[1])
        center_source = "execute"
    else:
        xc_auto, yc_auto = estimate_hpc_center(model)
        center_source = "lonC/latC"

    model_w_arcsec, model_h_arcsec = infer_fov_from_execute(
        loader_name=loader,
        model_path=model_path,
        dsun_cm=float(model["DSun"][0]),
        fallback_nx=int(model["Nx"][0]),
        fallback_ny=int(model["Ny"][0]),
        fallback_dx_cm=float(model["dx"][0]),
    )
    return model, model_dt, center_source, xc_auto, yc_auto, float(model_w_arcsec), float(model_h_arcsec), applied_overrides


def resolve_box_from_args(
    args: Any,
    *,
    xc_auto: float,
    yc_auto: float,
    model_w_arcsec: float,
    model_h_arcsec: float,
) -> tuple[float, float, float, float, int, int, float, float]:
    dx = float(args.dx) if args.dx is not None else None
    dy = float(args.dy) if args.dy is not None else None
    if args.pixel_scale_arcsec is not None:
        if dx is None:
            dx = float(args.pixel_scale_arcsec)
        if dy is None:
            dy = float(args.pixel_scale_arcsec)
    if dx is None:
        dx = 2.0
    if dy is None:
        dy = 2.0
    if dx <= 0 or dy <= 0:
        raise ValueError("--dx and --dy must be positive.")

    if args.xrange is not None:
        xmin, xmax = float(args.xrange[0]), float(args.xrange[1])
        if xmax <= xmin:
            raise ValueError("--xrange requires XMAX > XMIN.")
        if getattr(args, "nx", None) is not None:
            print("Note: --xrange provided, so --nx is ignored and recomputed from --dx.")
        xc = 0.5 * (xmin + xmax)
        nx = max(16, int(np.ceil((xmax - xmin) / dx)))
    else:
        xc = float(args.xc) if args.xc is not None else float(xc_auto)
        nx = max(16, int(args.nx)) if args.nx is not None else max(16, int(np.ceil(model_w_arcsec / dx)))

    if args.yrange is not None:
        ymin, ymax = float(args.yrange[0]), float(args.yrange[1])
        if ymax <= ymin:
            raise ValueError("--yrange requires YMAX > YMIN.")
        if getattr(args, "ny", None) is not None:
            print("Note: --yrange provided, so --ny is ignored and recomputed from --dy.")
        yc = 0.5 * (ymin + ymax)
        ny = max(16, int(np.ceil((ymax - ymin) / dy)))
    else:
        yc = float(args.yc) if args.yc is not None else float(yc_auto)
        ny = max(16, int(args.ny)) if args.ny is not None else max(16, int(np.ceil(model_h_arcsec / dy)))

    fov_x = float(nx) * float(dx)
    fov_y = float(ny) * float(dy)
    return float(xc), float(yc), float(dx), float(dy), int(nx), int(ny), fov_x, fov_y


def plasma_defaults() -> PlasmaDefaults:
    # Coronal plasma parameters aligned with IDL RenderExample* defaults.
    w = np.array([1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4], dtype=np.float64)
    shtable = np.outer(w, w)
    shtable[6, 6] = 0.1
    return PlasmaDefaults(
        tbase=1e6,
        nbase=1e8,
        q0=0.0217,
        a=0.3,
        b=2.7,
        shtable=shtable,
    )


def apply_model_observer_overrides(
    model: Any,
    *,
    dsun_cm: float | None = None,
    lonc_deg: float | None = None,
    b0sun_deg: float | None = None,
) -> dict[str, float]:
    applied: dict[str, float] = {}
    if dsun_cm is not None:
        model["DSun"] = float(dsun_cm)
        applied["DSun_cm"] = float(dsun_cm)
    if lonc_deg is not None:
        model["lonC"] = float(lonc_deg)
        applied["lonC_deg"] = float(lonc_deg)
    if b0sun_deg is not None:
        model["b0Sun"] = float(b0sun_deg)
        applied["b0Sun_deg"] = float(b0sun_deg)
    return applied


def model_obstime_iso(model: Any) -> str:
    return Time(float(model["obstime"][0]) + 283996800.0, format="unix").isot


def prepare_common_inputs(
    args: Any,
    *,
    prefer_execute_center: bool = True,
    observer_overrides: dict[str, float | None] | None = None,
) -> CommonRenderInputs:
    os.environ["OMP_NUM_THREADS"] = str(args.omp_threads)

    model_path = Path(args.model_path)
    loader = resolve_model_loader(model_path, args.model_format)
    ebtel_path, ebtel_c, ebtel_dt = load_ebtel_from_arg(args.ebtel_path)
    model, model_dt, center_source, xc_auto, yc_auto, model_w_arcsec, model_h_arcsec, applied_overrides = load_model_and_fov(
        model_path,
        loader,
        prefer_execute_center=prefer_execute_center,
        observer_overrides=observer_overrides,
    )
    xc, yc, dx, dy, nx, ny, fov_x, fov_y = resolve_box_from_args(
        args,
        xc_auto=xc_auto,
        yc_auto=yc_auto,
        model_w_arcsec=model_w_arcsec,
        model_h_arcsec=model_h_arcsec,
    )

    return CommonRenderInputs(
        model_path=model_path,
        loader=loader,
        model=model,
        model_dt=model_dt,
        ebtel_path=ebtel_path,
        ebtel_c=ebtel_c,
        ebtel_dt=ebtel_dt,
        center_source=center_source,
        xc=xc,
        yc=yc,
        dx=dx,
        dy=dy,
        nx=nx,
        ny=ny,
        fov_x=fov_x,
        fov_y=fov_y,
        observer_overrides_applied=applied_overrides,
    )
