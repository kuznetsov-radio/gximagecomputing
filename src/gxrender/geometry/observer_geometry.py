from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from pyampp.geometry import (
    GeometryContract,
    build_ephemeris_from_pb0r as pyampp_build_ephemeris_from_pb0r,
    build_fov_box_from_red_box_world,
    build_pb0r_metadata_from_ephemeris,
    compute_inscribing_fov_from_world,
    make_observer_wcs_header,
    normalize_observer_key,
    resolve_named_observer,
    resolve_observer_with_info,
    world_corners_from_geometry_contract,
)
from sunpy.coordinates import frames, get_horizons_coord
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst

_RSUN_METERS = 695700000.0
_HORIZONS_TARGETS = {
    "solar orbiter": "Solar Orbiter",
    "stereo-a": "STEREO-A",
    "stereo-b": "STEREO-B",
}


@dataclass(frozen=True)
class ResolvedObserverGeometry:
    observer_name: str
    l0_deg: float
    b0_deg: float
    dsun_cm: float
    render_lonc_deg: float
    render_b0_deg: float
    render_dsun_cm: float
    observer_source: str
    warnings: tuple[str, ...]
    rsun_cm: float | None = None
    rsun_arcsec: float | None = None


@dataclass(frozen=True)
class MetadataObserverState:
    observer_name: str
    l0_deg: float
    b0_deg: float
    dsun_cm: float
    source: str
    rsun_cm: float | None = None
    rsun_arcsec: float | None = None


def model_time_from_model(model: Any) -> Time:
    return Time(float(model["obstime"][0]) + 283996800.0, format="unix")


def normalize_observer_name(name: str | None) -> str | None:
    if name is None:
        return None
    normalized = normalize_observer_key(name)
    cleaned = str(normalized).strip().lower()
    return cleaned or None


def _pretty_observer_name(name: str | None) -> str:
    if not name:
        return "custom"
    parts = str(name).replace("-", " ").split()
    return " ".join(p.capitalize() for p in parts)


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        num = float(np.asarray(value).reshape(-1)[0])
    except Exception:
        try:
            num = float(value)
        except Exception:
            return None
    return num if np.isfinite(num) else None


def _normalize_dsun_cm(value: Any) -> float | None:
    dsun = _as_float(value)
    if dsun is None:
        return None
    return dsun * 100.0 if dsun < 1e12 else dsun


def _normalize_rsun_cm(value: Any) -> float | None:
    rsun = _as_float(value)
    if rsun is None:
        return None
    return rsun * 100.0 if rsun < 1e10 else rsun


def _normalize_rsun_arcsec(value: Any) -> float | None:
    rsun = _as_float(value)
    if rsun is None or not np.isfinite(rsun) or rsun <= 0:
        return None
    return float(rsun)


def _finalize_rsun_state(
    *,
    dsun_cm: float,
    rsun_cm: float | None = None,
    rsun_arcsec: float | None = None,
) -> tuple[float | None, float | None]:
    rsun_cm_value = _normalize_rsun_cm(rsun_cm)
    if rsun_cm_value is None:
        rsun_cm_value = _RSUN_METERS * 100.0
    rsun_arcsec_value = _normalize_rsun_arcsec(rsun_arcsec)
    if rsun_arcsec_value is None and np.isfinite(dsun_cm) and dsun_cm > 0:
        ratio = np.clip(float(rsun_cm_value) / float(dsun_cm), -1.0, 1.0)
        rsun_arcsec_value = float(np.arcsin(ratio) * u.rad.to(u.arcsec))
    return rsun_cm_value, rsun_arcsec_value


def _metadata_lookup(model_metadata: dict[str, Any] | None, *names: str) -> Any:
    if not isinstance(model_metadata, dict):
        return None
    lower_map = {str(k).lower(): v for k, v in model_metadata.items()}
    for name in names:
        key = str(name).lower()
        if key in lower_map:
            return lower_map[key]
    return None


def _nested_lookup(payload: dict[str, Any] | None, *path: str) -> Any:
    current: Any = payload
    for part in path:
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def _model_render_triad(model: Any) -> tuple[float, float, float]:
    return (
        float(model["lonC"][0]),
        float(model["b0Sun"][0]),
        float(model["DSun"][0]),
    )


def _render_lonc_for_observer(
    model_metadata: dict[str, Any] | None,
    *,
    observer_lon_deg: float,
    observer_lat_deg: float,
    observer_dsun_cm: float,
    model_time: Time,
    fallback_lonc_deg: float,
) -> float:
    lon_ref = _as_float(_metadata_lookup(model_metadata, "lon", "crval1"))
    if lon_ref is None:
        return float(fallback_lonc_deg)
    try:
        observer = SkyCoord(
            lon=float(observer_lon_deg) * u.deg,
            lat=float(observer_lat_deg) * u.deg,
            radius=float(observer_dsun_cm) * u.cm,
            frame=frames.HeliographicStonyhurst,
            obstime=model_time,
        )
        observer_hgc = observer.transform_to(frames.HeliographicCarrington(obstime=model_time, observer="self"))
        return float(lon_ref - observer_hgc.lon.to_value(u.deg))
    except Exception:
        return float(fallback_lonc_deg)


def _coerce_hgs(coord: SkyCoord, model_time: Time) -> tuple[float, float, float]:
    hgs = coord.transform_to(frames.HeliographicStonyhurst(obstime=model_time))
    return (
        float(hgs.lon.to_value(u.deg)),
        float(hgs.lat.to_value(u.deg)),
        float(hgs.radius.to_value(u.cm)),
    )


def _resolve_named_observer_hgs(name: str, model_time: Time) -> tuple[float, float, float] | None:
    coord = resolve_named_observer(name, model_time)
    if coord is not None:
        return _coerce_hgs(coord, model_time)

    # Compatibility fallback for environments where pyAMPP cannot resolve
    # remote ephemeris but direct SunPy paths are still available.
    try:
        coord = get_body_heliographic_stonyhurst(name, model_time)
    except Exception:
        target = _HORIZONS_TARGETS.get(name, name)
        try:
            coord = get_horizons_coord(target, model_time)
        except Exception:
            return None
    return _coerce_hgs(coord, model_time)


def _resolve_metadata_observer_hgs(
    model_time: Time,
    observer_metadata: dict[str, Any] | None,
    observer_name: str | None,
) -> tuple[str, tuple[float, float, float], str | None] | None:
    if not isinstance(observer_metadata, dict):
        return None

    has_observer_payload = any(key in observer_metadata for key in ("name", "ephemeris", "pb0r", "fov", "fov_box"))
    requested = observer_name or _nested_lookup(observer_metadata, "name")
    if not has_observer_payload and requested is None:
        return None

    context: dict[str, Any] = {"observer": dict(observer_metadata)}
    coord, warning, used_key = resolve_observer_with_info(context, requested, model_time)
    if coord is None:
        return None
    return str(used_key or normalize_observer_key(requested or "earth")), _coerce_hgs(coord, model_time), warning


def _metadata_square_fov(observer_metadata: dict[str, Any] | None) -> bool:
    for path in (("fov", "square"), ("fov_box", "square")):
        value = _nested_lookup(observer_metadata, *path)
        if value is not None:
            return bool(value)
    return False


def _resolve_saved_observer_from_model_metadata(model_metadata: dict[str, Any] | None) -> MetadataObserverState | None:
    name = normalize_observer_name(_metadata_lookup(model_metadata, "observer_name", "observer"))
    if not name:
        return None

    b0_deg = _as_float(_metadata_lookup(model_metadata, "observer_b0_deg"))
    l0_deg = _as_float(_metadata_lookup(model_metadata, "observer_l0_deg"))
    rsun_cm = _normalize_rsun_cm(_metadata_lookup(model_metadata, "observer_rsun_cm"))
    rsun_arcsec = _normalize_rsun_arcsec(_metadata_lookup(model_metadata, "observer_rsun_arcsec"))
    dsun_cm = _normalize_dsun_cm(_metadata_lookup(model_metadata, "observer_dsun_cm", "dsun_obs"))

    if dsun_cm is None and rsun_cm is not None and rsun_arcsec is not None:
        sin_rsun = np.sin(float(rsun_arcsec) * u.arcsec.to(u.rad))
        if np.isfinite(sin_rsun) and sin_rsun > 0:
            dsun_cm = float(rsun_cm / sin_rsun)

    if b0_deg is None or l0_deg is None or dsun_cm is None:
        return None

    return MetadataObserverState(
        observer_name=name,
        l0_deg=float(l0_deg),
        b0_deg=float(b0_deg),
        dsun_cm=float(dsun_cm),
        source="saved_observer_metadata",
        rsun_cm=rsun_cm,
        rsun_arcsec=rsun_arcsec,
    )


def _resolve_carrington_observer_from_model_metadata(
    model_time: Time,
    model_metadata: dict[str, Any] | None,
) -> MetadataObserverState | None:
    crln_obs = _as_float(_metadata_lookup(model_metadata, "crln_obs"))
    crlt_obs = _as_float(_metadata_lookup(model_metadata, "crlt_obs", "hglt_obs"))
    dsun_cm = _normalize_dsun_cm(_metadata_lookup(model_metadata, "dsun_obs", "observer_dsun_cm"))
    if crln_obs is None or crlt_obs is None or dsun_cm is None:
        return None

    try:
        observer_hgc = SkyCoord(
            lon=float(crln_obs) * u.deg,
            lat=float(crlt_obs) * u.deg,
            radius=float(dsun_cm) * u.cm,
            frame=frames.HeliographicCarrington(observer="self", obstime=model_time),
        )
        observer_hgs = observer_hgc.transform_to(frames.HeliographicStonyhurst(obstime=model_time))
    except Exception:
        return None

    return MetadataObserverState(
        observer_name="custom",
        l0_deg=float(observer_hgs.lon.to_value(u.deg)),
        b0_deg=float(observer_hgs.lat.to_value(u.deg)),
        dsun_cm=float(observer_hgs.radius.to_value(u.cm)),
        source="model_metadata_carrington",
    )


def should_use_saved_observer_fov(
    observer_metadata: dict[str, Any] | None,
    *,
    resolved_observer_name: str,
) -> bool:
    """Return whether saved FOV metadata should be used.

    Validates that the observer metadata contains FOV data. Does NOT perform
    observer name matching - that validation is delegated to pyAMPP's observer
    and FOV APIs which are the authoritative source.

    If pyAMPP loaded this observer metadata, then both the observer identity
    and FOV are already validated by pyAMPP.
    """
    if not isinstance(observer_metadata, dict):
        return False

    # pyAMPP observer blocks have "fov" key with the FOV metadata
    saved = observer_metadata.get("fov")
    if not isinstance(saved, dict):
        return False

    # If FOV metadata exists, we should use it. pyAMPP ensures consistency
    # between observer identity and FOV metadata.
    return True


def _cli_render_triad(cli_args: Any, model: Any) -> tuple[float, float, float]:
    model_lonc, model_b0, model_dsun = _model_render_triad(model)
    dsun = _normalize_dsun_cm(getattr(cli_args, "dsun_cm", None))
    if dsun is None:
        dsun = model_dsun
    lonc = _as_float(getattr(cli_args, "lonc_deg", None))
    if lonc is None:
        lonc = model_lonc
    b0sun = _as_float(getattr(cli_args, "b0sun_deg", None))
    if b0sun is None:
        b0sun = model_b0
    return float(lonc), float(b0sun), float(dsun)


def resolve_observer_geometry(
    model: Any,
    cli_args: Any,
    model_metadata: dict[str, Any] | None,
    observer_metadata: dict[str, Any] | None = None,
) -> ResolvedObserverGeometry:
    warnings: list[str] = []
    model_time = model_time_from_model(model)
    model_lonc, model_b0, model_dsun = _model_render_triad(model)

    cli_name = normalize_observer_name(getattr(cli_args, "observer", None))
    if cli_name:
        resolved = _resolve_named_observer_hgs(cli_name, model_time)
        if resolved is not None:
            l0_deg, b0_deg, dsun_cm = resolved
            render_lonc = _render_lonc_for_observer(
                model_metadata,
                observer_lon_deg=l0_deg,
                observer_lat_deg=b0_deg,
                observer_dsun_cm=dsun_cm,
                model_time=model_time,
                fallback_lonc_deg=model_lonc,
            )
            rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=dsun_cm)
            return ResolvedObserverGeometry(
                cli_name,
                l0_deg,
                b0_deg,
                dsun_cm,
                render_lonc,
                b0_deg,
                dsun_cm,
                "cli_observer",
                tuple(warnings),
                rsun_cm=rsun_cm,
                rsun_arcsec=rsun_arcsec,
            )
        warnings.append(f"CLI observer '{getattr(cli_args, 'observer', cli_name)}' could not be resolved")

    any_cli_triad = any(getattr(cli_args, key, None) is not None for key in ("lonc_deg", "b0sun_deg", "dsun_cm"))
    cli_lonc, cli_b0, cli_dsun = _cli_render_triad(cli_args, model)
    if any_cli_triad:
        actual_name = "earth"
        rsun_cm = None
        rsun_arcsec = None
        resolved_meta = _resolve_metadata_observer_hgs(model_time, observer_metadata, None)
        if resolved_meta is not None:
            actual_name, (l0_deg, b0_deg, dsun_cm), warning = resolved_meta
            if warning:
                warnings.append(str(warning))
        else:
            resolved_earth = _resolve_named_observer_hgs("earth", model_time)
            if resolved_earth is None:
                raise ValueError("Could not resolve Earth observer from pyAMPP geometry API")
            l0_deg, b0_deg, dsun_cm = resolved_earth
        rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=dsun_cm, rsun_cm=rsun_cm, rsun_arcsec=rsun_arcsec)
        return ResolvedObserverGeometry(
            actual_name,
            l0_deg,
            b0_deg,
            dsun_cm,
            cli_lonc,
            cli_b0,
            cli_dsun,
            "cli_triad",
            tuple(warnings),
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    resolved_meta = _resolve_metadata_observer_hgs(model_time, observer_metadata, None)
    if resolved_meta is not None:
        observer_name, (l0_deg, b0_deg, dsun_cm), warning = resolved_meta
        if warning:
            warnings.append(str(warning))
        render_lonc = _render_lonc_for_observer(
            model_metadata,
            observer_lon_deg=l0_deg,
            observer_lat_deg=b0_deg,
            observer_dsun_cm=dsun_cm,
            model_time=model_time,
            fallback_lonc_deg=model_lonc,
        )
        rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=dsun_cm)
        return ResolvedObserverGeometry(
            observer_name,
            l0_deg,
            b0_deg,
            dsun_cm,
            render_lonc,
            b0_deg,
            dsun_cm,
            "saved_observer_metadata",
            tuple(warnings),
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    saved_state = _resolve_saved_observer_from_model_metadata(model_metadata)
    if saved_state is not None:
        render_lonc = _render_lonc_for_observer(
            model_metadata,
            observer_lon_deg=saved_state.l0_deg,
            observer_lat_deg=saved_state.b0_deg,
            observer_dsun_cm=saved_state.dsun_cm,
            model_time=model_time,
            fallback_lonc_deg=model_lonc,
        )
        rsun_cm, rsun_arcsec = _finalize_rsun_state(
            dsun_cm=saved_state.dsun_cm,
            rsun_cm=saved_state.rsun_cm,
            rsun_arcsec=saved_state.rsun_arcsec,
        )
        return ResolvedObserverGeometry(
            saved_state.observer_name,
            saved_state.l0_deg,
            saved_state.b0_deg,
            saved_state.dsun_cm,
            render_lonc,
            saved_state.b0_deg,
            saved_state.dsun_cm,
            saved_state.source,
            tuple(warnings),
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    carrington_state = _resolve_carrington_observer_from_model_metadata(model_time, model_metadata)
    if carrington_state is not None:
        render_lonc = _render_lonc_for_observer(
            model_metadata,
            observer_lon_deg=carrington_state.l0_deg,
            observer_lat_deg=carrington_state.b0_deg,
            observer_dsun_cm=carrington_state.dsun_cm,
            model_time=model_time,
            fallback_lonc_deg=model_lonc,
        )
        rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=carrington_state.dsun_cm)
        return ResolvedObserverGeometry(
            carrington_state.observer_name,
            carrington_state.l0_deg,
            carrington_state.b0_deg,
            carrington_state.dsun_cm,
            render_lonc,
            carrington_state.b0_deg,
            carrington_state.dsun_cm,
            carrington_state.source,
            tuple(warnings),
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    meta_name = normalize_observer_name(
        _nested_lookup(observer_metadata, "name")
        or _metadata_lookup(model_metadata, "observer_name", "observer", "observatory", "obsrvtry")
    )
    if meta_name:
        resolved = _resolve_named_observer_hgs(meta_name, model_time)
        if resolved is not None:
            l0_deg, b0_deg, dsun_cm = resolved
            render_lonc = _render_lonc_for_observer(
                model_metadata,
                observer_lon_deg=l0_deg,
                observer_lat_deg=b0_deg,
                observer_dsun_cm=dsun_cm,
                model_time=model_time,
                fallback_lonc_deg=model_lonc,
            )
            rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=dsun_cm)
            return ResolvedObserverGeometry(
                meta_name,
                l0_deg,
                b0_deg,
                dsun_cm,
                render_lonc,
                b0_deg,
                dsun_cm,
                "model_metadata_observer",
                tuple(warnings),
                rsun_cm=rsun_cm,
                rsun_arcsec=rsun_arcsec,
            )
        warnings.append(f"Model observer '{meta_name}' could not be resolved")

    resolved_earth = _resolve_named_observer_hgs("earth", model_time)
    if resolved_earth is None:
        raise ValueError("Could not resolve default Earth observer from pyAMPP geometry API")
    l0_deg, b0_deg, dsun_cm = resolved_earth
    rsun_cm, rsun_arcsec = _finalize_rsun_state(dsun_cm=dsun_cm)
    return ResolvedObserverGeometry(
        "earth",
        l0_deg,
        b0_deg,
        dsun_cm,
        model_lonc,
        model_b0,
        model_dsun,
        "default_earth",
        tuple(warnings),
        rsun_cm=rsun_cm,
        rsun_arcsec=rsun_arcsec,
    )


def observer_summary(geometry: ResolvedObserverGeometry) -> str:
    if geometry.observer_name != "custom":
        return _pretty_observer_name(geometry.observer_name)
    return f"custom (l0={geometry.l0_deg:.3f}, b0={geometry.b0_deg:.3f}, dsun={geometry.dsun_cm:.6g} cm)"


def build_observer_coordinate(geometry: ResolvedObserverGeometry, obs_time: str | Time) -> SkyCoord:
    obstime = Time(obs_time)
    return SkyCoord(
        lon=float(geometry.l0_deg) * u.deg,
        lat=float(geometry.b0_deg) * u.deg,
        radius=float(geometry.dsun_cm) * u.cm,
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
    )


def _geometry_contract_from_metadata(
    model: Any,
    model_metadata: dict[str, Any] | None,
    *,
    obstime: Time,
) -> GeometryContract | None:
    if not isinstance(model_metadata, dict):
        model_metadata = {}

    value = model_metadata.get("geometry_contract")
    if isinstance(value, GeometryContract):
        return value
    if isinstance(value, dict):
        try:
            return GeometryContract.from_dict(value)
        except Exception:
            pass

    try:
        nx = int(_as_float(model_metadata.get("box_nx")) or int(model["Nx"][0]))
        ny = int(_as_float(model_metadata.get("box_ny")) or int(model["Ny"][0]))
        dz = np.asarray(model["dz"][0], dtype=np.float64)
        nz = int(_as_float(model_metadata.get("box_nz")) or dz.shape[0])
        rsun_cm = float(model["RSun"][0])

        dr_x = float(_as_float(model_metadata.get("box_dr_x")) or (float(model["dx"][0]) / rsun_cm))
        dr_y = float(_as_float(model_metadata.get("box_dr_y")) or (float(model["dy"][0]) / rsun_cm))
        if _as_float(model_metadata.get("box_dr_z")) is not None:
            dr_z = float(_as_float(model_metadata.get("box_dr_z")))
        else:
            dz_med = float(np.nanmedian(dz[np.isfinite(dz)])) if np.any(np.isfinite(dz)) else float(model["dx"][0])
            dr_z = float(dz_med / rsun_cm)

        rsun_ref_m = _as_float(model_metadata.get("rsun_ref"))
        if rsun_ref_m is None:
            rsun_ref_m = rsun_cm / 100.0
        anchor_lon = _as_float(_metadata_lookup(model_metadata, "lon", "crval1"))
        if anchor_lon is None:
            anchor_lon = float(model["lonC"][0])
        anchor_lat = _as_float(_metadata_lookup(model_metadata, "lat", "crval2"))
        if anchor_lat is None:
            anchor_lat = float(model["latC"][0]) if "latC" in model.dtype.names else 0.0
        return GeometryContract(
            nx=nx,
            ny=ny,
            nz=nz,
            dr_x=dr_x,
            dr_y=dr_y,
            dr_z=dr_z,
            rsun_m=float(rsun_ref_m),
            anchor_lon_deg=float(anchor_lon),
            anchor_lat_deg=float(anchor_lat),
            anchor_radius_rsun=1.0,
            frame="heliographic_stonyhurst",
            obstime=obstime.isot,
            inferred_from="gxrender_metadata_fallback",
        )
    except Exception:
        return None


def _world_corners_from_model_metadata(
    model: Any,
    model_metadata: dict[str, Any] | None,
    obstime: Time,
    observer: str | SkyCoord | None,
) -> SkyCoord:
    contract = _geometry_contract_from_metadata(model, model_metadata, obstime=obstime)
    if contract is None:
        raise ValueError("Missing geometry_contract metadata; cannot delegate observer/FOV geometry to pyAMPP")
    world = world_corners_from_geometry_contract(contract, obstime=obstime, observer=observer)
    if world is None:
        raise ValueError("pyAMPP could not build world corners from geometry_contract")
    return world


def compute_projected_fov_for_observer(
    model: Any,
    observer_geometry: ResolvedObserverGeometry,
    *,
    model_metadata: dict[str, Any] | None = None,
    observer_metadata: dict[str, Any] | None = None,
) -> tuple[float, float, float, float]:
    fov = compute_inscribing_fov(
        model,
        observer_geometry,
        model_metadata=model_metadata,
        observer_metadata=observer_metadata,
    )
    return (
        float(fov["xc_arcsec"]),
        float(fov["yc_arcsec"]),
        float(fov["xsize_arcsec"]),
        float(fov["ysize_arcsec"]),
    )


def compute_inscribing_fov(
    model: Any,
    observer_geometry: ResolvedObserverGeometry,
    *,
    model_metadata: dict[str, Any] | None = None,
    observer_metadata: dict[str, Any] | None = None,
    pad_arcsec: float = 0.0,
) -> dict[str, Any]:
    obstime = model_time_from_model(model)
    observer = build_observer_coordinate(observer_geometry, obstime)
    corners_hgs = _world_corners_from_model_metadata(model, model_metadata, obstime, observer)
    fov = compute_inscribing_fov_from_world(
        corners_hgs,
        observer=observer,
        obstime=obstime,
        pad_arcsec=max(float(pad_arcsec), 0.0),
    )
    if fov is None:
        raise ValueError("pyAMPP failed to compute inscribing observer FOV")

    xmin = float(fov["xmin_arcsec"])
    xmax = float(fov["xmax_arcsec"])
    ymin = float(fov["ymin_arcsec"])
    ymax = float(fov["ymax_arcsec"])
    if _metadata_square_fov(observer_metadata):
        side = max(xmax - xmin, ymax - ymin)
        xc = 0.5 * (xmin + xmax)
        yc = 0.5 * (ymin + ymax)
        xmin = xc - 0.5 * side
        xmax = xc + 0.5 * side
        ymin = yc - 0.5 * side
        ymax = yc + 0.5 * side

    return {
        "xc_arcsec": 0.5 * (xmin + xmax),
        "yc_arcsec": 0.5 * (ymin + ymax),
        "xsize_arcsec": xmax - xmin,
        "ysize_arcsec": ymax - ymin,
        "xmin_arcsec": xmin,
        "xmax_arcsec": xmax,
        "ymin_arcsec": ymin,
        "ymax_arcsec": ymax,
        "corners_hpc": fov.get("corners_hpc"),
    }


def compute_inscribing_fov_box(
    model: Any,
    observer_geometry: ResolvedObserverGeometry,
    *,
    model_metadata: dict[str, Any] | None = None,
    observer_metadata: dict[str, Any] | None = None,
    pad_xy_arcsec: float = 0.0,
    pad_z_frac: float = 0.10,
) -> dict[str, Any]:
    obstime = model_time_from_model(model)
    observer = build_observer_coordinate(observer_geometry, obstime)
    corners_hgs = _world_corners_from_model_metadata(model, model_metadata, obstime, observer)
    fov_box = build_fov_box_from_red_box_world(
        corners_hgs,
        observer=observer,
        obstime=obstime,
        pad_xy_arcsec=max(float(pad_xy_arcsec), 0.0),
        pad_z_frac=max(float(pad_z_frac), 0.0),
    )
    if fov_box is None:
        raise ValueError("pyAMPP failed to compute inscribing observer FOV box")

    if _metadata_square_fov(observer_metadata):
        side = max(float(fov_box["xsize_arcsec"]), float(fov_box["ysize_arcsec"]))
        fov_box = dict(fov_box)
        fov_box["xsize_arcsec"] = side
        fov_box["ysize_arcsec"] = side
        fov_box["xmin_arcsec"] = float(fov_box["xc_arcsec"]) - 0.5 * side
        fov_box["xmax_arcsec"] = float(fov_box["xc_arcsec"]) + 0.5 * side
        fov_box["ymin_arcsec"] = float(fov_box["yc_arcsec"]) - 0.5 * side
        fov_box["ymax_arcsec"] = float(fov_box["yc_arcsec"]) + 0.5 * side
    return dict(fov_box)


def resolve_simbox_from_observer_and_model(
    *,
    explicit_xc: float | None = None,
    explicit_yc: float | None = None,
    explicit_xsize: float | None = None,
    explicit_ysize: float | None = None,
    saved_fov: dict[str, Any] | None = None,
    computed_fov: dict[str, Any] | None = None,
) -> tuple[str, float, float, float, float] | None:
    if None not in (explicit_xc, explicit_yc, explicit_xsize, explicit_ysize):
        return ("explicit", float(explicit_xc), float(explicit_yc), float(explicit_xsize), float(explicit_ysize))
    if isinstance(saved_fov, dict):
        xc = _as_float(saved_fov.get("xc_arcsec"))
        yc = _as_float(saved_fov.get("yc_arcsec"))
        xsize = _as_float(saved_fov.get("xsize_arcsec"))
        ysize = _as_float(saved_fov.get("ysize_arcsec"))
        if None not in (xc, yc, xsize, ysize):
            return ("saved_observer_fov", float(xc), float(yc), float(xsize), float(ysize))
    if isinstance(computed_fov, dict):
        xc = _as_float(computed_fov.get("xc_arcsec"))
        yc = _as_float(computed_fov.get("yc_arcsec"))
        xsize = _as_float(computed_fov.get("xsize_arcsec"))
        ysize = _as_float(computed_fov.get("ysize_arcsec"))
        if None not in (xc, yc, xsize, ysize):
            return ("inscribing_fov", float(xc), float(yc), float(xsize), float(ysize))
    return None


def compute_sunpy_wcs_header(
    *,
    nx: int,
    ny: int,
    xc_arcsec: float,
    yc_arcsec: float,
    dx_arcsec: float,
    dy_arcsec: float,
    obs_time: str | Time,
    observer_geometry: ResolvedObserverGeometry,
    bunit: str,
) -> fits.Header:
    observer = build_observer_coordinate(observer_geometry, obs_time)
    rsun_ref_m = (
        float(observer_geometry.rsun_cm) / 100.0
        if observer_geometry.rsun_cm is not None and np.isfinite(observer_geometry.rsun_cm)
        else _RSUN_METERS
    )
    return make_observer_wcs_header(
        nx=int(nx),
        ny=int(ny),
        xc_arcsec=float(xc_arcsec),
        yc_arcsec=float(yc_arcsec),
        dx_arcsec=float(dx_arcsec),
        dy_arcsec=float(dy_arcsec),
        observer=observer,
        obs_time=Time(obs_time).isot,
        bunit=str(bunit),
        observer_name=observer_geometry.observer_name,
        rsun_ref_m=rsun_ref_m,
        rsun_obs_arcsec=observer_geometry.rsun_arcsec,
    )


def build_ephemeris_from_pb0r(
    *,
    b0_deg: Any,
    l0_deg: Any,
    rsun_arcsec: Any,
    obs_date: str | Time | None = None,
    rsun_cm: Any = None,
) -> dict[str, Any] | None:
    return pyampp_build_ephemeris_from_pb0r(
        b0_deg=b0_deg,
        l0_deg=l0_deg,
        rsun_arcsec=rsun_arcsec,
        obs_date=Time(obs_date).isot if obs_date is not None else None,
        rsun_cm=rsun_cm,
    )


def build_pb0r_from_ephemeris(ephemeris: dict[str, Any] | None) -> dict[str, Any] | None:
    return build_pb0r_metadata_from_ephemeris(ephemeris)
