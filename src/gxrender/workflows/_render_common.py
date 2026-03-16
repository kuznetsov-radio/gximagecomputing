from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.time import Time

from gxrender.io.ebtel import load_ebtel, load_ebtel_none
from gxrender.geometry.observer_geometry import (
    ResolvedObserverGeometry,
    compute_inscribing_fov,
    resolve_observer_geometry,
    resolve_simbox_from_observer_and_model,
)
from gxrender.io.model import (
    estimate_hpc_center,
    infer_center_from_execute,
    infer_fov_from_execute,
    load_model_hdf_with_observer,
    load_model_sav_with_observer,
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
    model_metadata: dict[str, Any]
    observer_geometry: ResolvedObserverGeometry
    observer_overrides_applied: dict[str, float]


@dataclass
class CoronalPlasmaParameters:
    tbase: float
    nbase: float
    q0: float
    a: float
    b: float
    mode: int
    selective_heating: bool
    shtable: np.ndarray | None


def _default_shtable() -> np.ndarray:
    weights = np.array([1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4], dtype=np.float64)
    shtable = np.outer(weights, weights)
    shtable[6, 6] = 0.1
    return shtable


def _coerce_shtable(value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if arr.shape != (7, 7):
        raise ValueError(f"SHtable must have shape (7, 7), got {arr.shape}.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("SHtable must contain only finite numeric values.")
    return arr


def _load_shtable_from_path(path_value: str | os.PathLike[str]) -> np.ndarray:
    path = Path(path_value).expanduser()
    suffix = path.suffix.lower()
    if suffix == ".npy":
        return _coerce_shtable(np.load(path, allow_pickle=False))
    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as npz:
            if len(npz.files) != 1:
                raise ValueError(f"SHtable NPZ must contain exactly one array, found {npz.files}.")
            return _coerce_shtable(npz[npz.files[0]])
    delimiter = "," if suffix == ".csv" else None
    return _coerce_shtable(np.loadtxt(path, delimiter=delimiter))


def resolve_model_loader(model_path: Path, model_format: str) -> str:
    if model_format == "auto":
        return "h5" if model_path.suffix.lower() in {".h5", ".hdf5"} else "sav"
    return str(model_format)


def load_ebtel_from_arg(ebtel_path_arg: str | None) -> tuple[str, Any, Any]:
    if ebtel_path_arg is None:
        raise ValueError(
            "EBTEL input must be explicit in the shared workflow. "
            'Provide a .sav path, or pass "" to disable DEM/DDM tables.'
        )

    ebtel_raw = str(ebtel_path_arg).strip()
    if ebtel_raw in {"", '""', "''"}:
        ebtel_path = ""
        ebtel_c, ebtel_dt = load_ebtel_none()
        return ebtel_path, ebtel_c, ebtel_dt

    ebtel_path = str(Path(ebtel_raw).expanduser())
    if not Path(ebtel_path).exists():
        raise FileNotFoundError(f"Explicit EBTEL file not found: {ebtel_path}")
    ebtel_c, ebtel_dt = load_ebtel(ebtel_path)
    return ebtel_path, ebtel_c, ebtel_dt


def load_model_and_fov(
    model_path: Path,
    loader: str,
    cli_args: Any,
    *,
    prefer_execute_center: bool = True,
) -> tuple[Any, Any, dict[str, Any], ResolvedObserverGeometry, str, float, float, float, float, dict[str, float]]:
    if loader == "h5":
        model, model_dt, model_metadata, observer_metadata = load_model_hdf_with_observer(str(model_path))
    else:
        model, model_dt, model_metadata, observer_metadata = load_model_sav_with_observer(str(model_path))

    observer_geometry = resolve_observer_geometry(model, cli_args, model_metadata, observer_metadata)
    applied_overrides = apply_model_observer_overrides(
        model,
        dsun_cm=observer_geometry.render_dsun_cm,
        lonc_deg=observer_geometry.render_lonc_deg,
        b0sun_deg=observer_geometry.render_b0_deg,
    )

    has_cli_observer_override = any(
        getattr(cli_args, key, None) is not None for key in ("observer", "lonc_deg", "b0sun_deg", "dsun_cm")
    )
    has_cli_view_override = any(
        getattr(cli_args, key, None) is not None
        for key in ("xc", "yc", "xrange", "yrange", "nx", "ny", "dx", "dy", "pixel_scale_arcsec")
    )

    prefer_saved_fov = not has_cli_observer_override and not has_cli_view_override
    if bool(getattr(cli_args, "auto_fov", False)):
        prefer_saved_fov = False
    if bool(getattr(cli_args, "use_saved_fov", False)):
        prefer_saved_fov = True

    saved_fov = None
    if isinstance(observer_metadata, dict) and prefer_saved_fov:
        saved_fov = observer_metadata.get("fov")
    computed_fov = None
    try:
        computed_fov = compute_inscribing_fov(
            model,
            observer_geometry,
            model_metadata=model_metadata,
            observer_metadata=observer_metadata,
        )
    except Exception:
        computed_fov = None

    simbox = resolve_simbox_from_observer_and_model(saved_fov=saved_fov, computed_fov=computed_fov)
    if simbox is not None:
        center_source, xc_auto, yc_auto, model_w_arcsec, model_h_arcsec = simbox
    else:
        center_exec = infer_center_from_execute(loader_name=loader, model_path=model_path) if prefer_execute_center else None
        if center_exec is not None and observer_geometry.observer_source == "default_earth":
            xc_auto, yc_auto = float(center_exec[0]), float(center_exec[1])
            center_source = "execute"
        else:
            xc_auto, yc_auto = estimate_hpc_center(
                model,
                observer_lon_deg=observer_geometry.l0_deg,
                observer_lat_deg=observer_geometry.b0_deg,
                observer_dsun_cm=observer_geometry.dsun_cm,
            )
            center_source = "lonC/latC"

        model_w_arcsec, model_h_arcsec = infer_fov_from_execute(
            loader_name=loader,
            model_path=model_path,
            dsun_cm=float(observer_geometry.render_dsun_cm),
            fallback_nx=int(model["Nx"][0]),
            fallback_ny=int(model["Ny"][0]),
            fallback_dx_cm=float(model["dx"][0]),
        )
    return (
        model,
        model_dt,
        dict(model_metadata),
        observer_geometry,
        center_source,
        xc_auto,
        yc_auto,
        float(model_w_arcsec),
        float(model_h_arcsec),
        applied_overrides,
    )


def _metadata_float(model_metadata: dict[str, Any], key: str) -> float | None:
    value = model_metadata.get(key)
    if value is None:
        return None
    try:
        scalar = float(np.asarray(value).reshape(-1)[0])
    except Exception:
        try:
            scalar = float(value)
        except Exception:
            return None
    return scalar if np.isfinite(scalar) else None


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
    if dx is None or dy is None:
        raise ValueError(
            "Output pixel scale must be explicit in the shared workflow. "
            "Provide --dx/--dy or --pixel-scale-arcsec."
        )
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


def resolve_plasma_parameters(args: Any) -> CoronalPlasmaParameters:
    missing = [
        name
        for name in ("tbase", "nbase", "q0", "a", "b")
        if getattr(args, name, None) is None
    ]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(
            "Coronal plasma inputs must be explicit in the shared workflow. "
            f"Missing: {joined}."
        )

    tbase = float(args.tbase)
    nbase = float(args.nbase)
    q0 = float(args.q0)
    a = float(args.a)
    b = float(args.b)
    mode_value = getattr(args, "corona_mode", None)
    mode = int(mode_value) if mode_value is not None else 0
    if getattr(args, "force_isothermal", False):
        mode |= 1
    if getattr(args, "interpol_b", False):
        mode |= 2
    if getattr(args, "analytical_nt", False):
        mode |= 4
    if mode < 0:
        raise ValueError("corona mode must be non-negative.")
    if tbase <= 0.0:
        raise ValueError("Tbase must be positive.")
    if nbase <= 0.0:
        raise ValueError("nbase must be positive.")
    if q0 < 0.0:
        raise ValueError("Q0 must be non-negative.")

    raw_shtable = getattr(args, "shtable", None)
    shtable_path = getattr(args, "shtable_path", None)
    selective_heating = bool(getattr(args, "selective_heating", False))
    if raw_shtable is not None:
        shtable = _coerce_shtable(raw_shtable)
        selective_heating = True
    elif shtable_path not in (None, ""):
        shtable = _load_shtable_from_path(shtable_path)
        selective_heating = True
    elif selective_heating:
        shtable = _default_shtable()
        setattr(args, "_default_shtable_applied", True)
        warnings.warn(
            "Selective heating was enabled without an explicit SHtable, so the standard 7x7 table was assumed.",
            stacklevel=2,
        )
    else:
        shtable = None

    return CoronalPlasmaParameters(
        tbase=tbase,
        nbase=nbase,
        q0=q0,
        a=a,
        b=b,
        mode=mode,
        selective_heating=selective_heating,
        shtable=shtable,
    )


def resolve_mw_frequencies(args: Any) -> list[float]:
    raw = getattr(args, "frequencies_ghz", None)
    if raw is None:
        raise ValueError(
            "Microwave frequency list must be explicit in the shared workflow. "
            "Provide --frequencies-ghz."
        )

    freqs = [float(v) for v in raw]
    if len(freqs) == 0:
        raise ValueError("At least one microwave frequency must be provided.")
    if any(v <= 0.0 for v in freqs):
        raise ValueError("Microwave frequencies must be positive.")
    return freqs


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
) -> CommonRenderInputs:
    os.environ["OMP_NUM_THREADS"] = str(args.omp_threads)

    model_path = Path(args.model_path)
    loader = resolve_model_loader(model_path, args.model_format)
    ebtel_path, ebtel_c, ebtel_dt = load_ebtel_from_arg(args.ebtel_path)
    (
        model,
        model_dt,
        model_metadata,
        observer_geometry,
        center_source,
        xc_auto,
        yc_auto,
        model_w_arcsec,
        model_h_arcsec,
        applied_overrides,
    ) = load_model_and_fov(
        model_path,
        loader,
        args,
        prefer_execute_center=prefer_execute_center,
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
        model_metadata=model_metadata,
        observer_geometry=observer_geometry,
        observer_overrides_applied=applied_overrides,
    )
