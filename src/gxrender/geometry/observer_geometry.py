from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from sunpy.coordinates import frames, get_horizons_coord
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst
from sunpy.map.header_helper import make_fitswcs_header
try:
    from sunpy.coordinates.screens import SphericalScreen
except Exception:  # pragma: no cover
    SphericalScreen = None

from gxrender.io.model import extract_geometry_from_execute

_RSUN_METERS = 695700000.0

_OBSERVER_ALIASES = {
    "earth": "earth",
    "terra": "earth",
    "solo": "solar orbiter",
    "solar orbiter": "solar orbiter",
    "solar-orbiter": "solar orbiter",
    "solarorbiter": "solar orbiter",
    "stereo a": "stereo-a",
    "stereo-a": "stereo-a",
    "stereoa": "stereo-a",
    "stereo ahead": "stereo-a",
    "stereo b": "stereo-b",
    "stereo-b": "stereo-b",
    "stereob": "stereo-b",
    "stereo behind": "stereo-b",
}

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
    if rsun_arcsec_value is None and dsun_cm is not None and np.isfinite(dsun_cm) and dsun_cm > 0:
        ratio = np.clip(float(rsun_cm_value) / float(dsun_cm), -1.0, 1.0)
        rsun_arcsec_value = float(np.arcsin(ratio) * u.rad.to(u.arcsec))
    return rsun_cm_value, rsun_arcsec_value


def model_time_from_model(model: Any) -> Time:
    return Time(float(model["obstime"][0]) + 283996800.0, format="unix")


def normalize_observer_name(name: str | None) -> str | None:
    if name is None:
        return None
    cleaned = " ".join(str(name).strip().split())
    if not cleaned:
        return None
    lowered = cleaned.lower()
    return _OBSERVER_ALIASES.get(lowered, lowered)


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


def _derive_dsun_cm_from_rsun_arcsec(rsun_arcsec: Any, rsun_cm: Any = None) -> float | None:
    rsun_arcsec_value = _normalize_rsun_arcsec(rsun_arcsec)
    if rsun_arcsec_value is None:
        return None
    rsun_cm_value = _normalize_rsun_cm(rsun_cm)
    if rsun_cm_value is None:
        rsun_cm_value = _RSUN_METERS * 100.0
    try:
        rsun_rad = float((rsun_arcsec_value * u.arcsec).to_value(u.rad))
        if not np.isfinite(rsun_rad) or rsun_rad <= 0:
            return None
        return float(rsun_cm_value / np.sin(rsun_rad))
    except Exception:
        return None


def build_ephemeris_from_pb0r(
    *,
    b0_deg: Any,
    l0_deg: Any,
    rsun_arcsec: Any,
    obs_date: str | Time | None = None,
    rsun_cm: Any = None,
) -> dict[str, Any] | None:
    b0_value = _as_float(b0_deg)
    l0_value = _as_float(l0_deg)
    rsun_arcsec_value = _normalize_rsun_arcsec(rsun_arcsec)
    if None in (b0_value, l0_value, rsun_arcsec_value):
        return None
    rsun_cm_value = _normalize_rsun_cm(rsun_cm)
    if rsun_cm_value is None:
        rsun_cm_value = _RSUN_METERS * 100.0
    dsun_cm = _derive_dsun_cm_from_rsun_arcsec(rsun_arcsec_value, rsun_cm_value)
    if dsun_cm is None:
        return None
    return {
        "obs_date": Time(obs_date).isot if obs_date is not None else None,
        "hgln_obs_deg": float(l0_value),
        "hglt_obs_deg": float(b0_value),
        "dsun_cm": float(dsun_cm),
        "rsun_cm": float(rsun_cm_value),
    }


def build_pb0r_from_ephemeris(ephemeris: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(ephemeris, dict):
        return None
    l0_deg = _as_float(ephemeris.get("hgln_obs_deg"))
    b0_deg = _as_float(ephemeris.get("hglt_obs_deg"))
    dsun_cm = _normalize_dsun_cm(ephemeris.get("dsun_cm"))
    rsun_cm = _normalize_rsun_cm(ephemeris.get("rsun_cm"))
    if None in (l0_deg, b0_deg, dsun_cm):
        return None
    rsun_arcsec = _normalize_rsun_arcsec(ephemeris.get("rsun_arcsec"))
    if rsun_arcsec is None and rsun_cm is not None:
        ratio = np.clip(float(rsun_cm) / float(dsun_cm), -1.0, 1.0)
        rsun_arcsec = float(np.arcsin(ratio) * u.rad.to(u.arcsec))
    return {
        "obs_date": ephemeris.get("obs_date"),
        "b0_deg": float(b0_deg),
        "l0_deg": float(l0_deg),
        "rsun_arcsec": rsun_arcsec,
    }


def _observer_from_sunpy(name: str, model_time: Time) -> tuple[float, float, float]:
    try:
        obs = get_body_heliographic_stonyhurst(name, model_time)
    except Exception:
        target = _HORIZONS_TARGETS.get(name, name)
        obs = get_horizons_coord(target, model_time)
    return (
        float(obs.lon.to_value(u.deg)),
        float(obs.lat.to_value(u.deg)),
        float(obs.radius.to_value(u.cm)),
    )


def _model_render_triad(model: Any) -> tuple[float, float, float]:
    return (
        float(model["lonC"][0]),
        float(model["b0Sun"][0]),
        float(model["DSun"][0]),
    )


def _carrington_observer_to_stonyhurst(
    carr_lon_deg: float,
    carr_lat_deg: float,
    dsun_cm: float,
    model_time: Time,
) -> tuple[float, float] | None:
    try:
        coord = SkyCoord(
            lon=float(carr_lon_deg) * u.deg,
            lat=float(carr_lat_deg) * u.deg,
            radius=float(dsun_cm) * u.cm,
            frame=frames.HeliographicCarrington,
            obstime=model_time,
            observer="self",
        )
        hgs = coord.transform_to(frames.HeliographicStonyhurst(obstime=model_time))
        return float(hgs.lon.to_value(u.deg)), float(hgs.lat.to_value(u.deg))
    except Exception:
        return None


def _nested_lookup(payload: dict[str, Any] | None, *path: str) -> Any:
    current: Any = payload
    for part in path:
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def _metadata_lookup(model_metadata: dict[str, Any] | None, *names: str) -> Any:
    if not isinstance(model_metadata, dict):
        return None
    lower_map = {str(k).lower(): v for k, v in model_metadata.items()}
    for name in names:
        key = str(name).lower()
        if key in lower_map:
            return lower_map[key]
    return None


def _metadata_observer_state(
    model_metadata: dict[str, Any] | None,
    observer_metadata: dict[str, Any] | None,
    model_time: Time,
) -> MetadataObserverState | None:
    observer_name = normalize_observer_name(
        _nested_lookup(observer_metadata, "name")
        or _metadata_lookup(model_metadata, "observer_name", "observer", "observatory", "obsrvtry")
    ) or "custom"

    ephemeris = _nested_lookup(observer_metadata, "ephemeris")
    pb0r = _nested_lookup(observer_metadata, "pb0r")

    rsun_cm = _normalize_rsun_cm(
        (ephemeris or {}).get("rsun_cm")
        if isinstance(ephemeris, dict)
        else _metadata_lookup(model_metadata, "observer_rsun_cm", "rsun_cm", "rsun_ref_cm", "rsun_ref_m")
    )
    rsun_arcsec = _normalize_rsun_arcsec(
        (pb0r or {}).get("rsun_arcsec")
        if isinstance(pb0r, dict)
        else _metadata_lookup(model_metadata, "observer_rsun_arcsec", "rsun_arcsec", "rsun_obs")
    )
    dsun_cm = _normalize_dsun_cm(
        (ephemeris or {}).get("dsun_cm")
        if isinstance(ephemeris, dict)
        else _metadata_lookup(model_metadata, "observer_dsun_cm", "dsun_obs", "dsun_cm", "dsun_m")
    )
    if dsun_cm is None:
        dsun_cm = _derive_dsun_cm_from_rsun_arcsec(rsun_arcsec, rsun_cm)

    l0_deg = _as_float(
        (pb0r or {}).get("l0_deg")
        if isinstance(pb0r, dict)
        else _metadata_lookup(model_metadata, "observer_hgln_obs_deg", "observer_l0_deg", "hgln_obs")
    )
    b0_deg = _as_float(
        (pb0r or {}).get("b0_deg")
        if isinstance(pb0r, dict)
        else _metadata_lookup(model_metadata, "observer_hglt_obs_deg", "observer_b0_deg", "hglt_obs")
    )
    if l0_deg is not None and b0_deg is not None and dsun_cm is not None:
        return MetadataObserverState(
            observer_name=observer_name,
            l0_deg=float(l0_deg),
            b0_deg=float(b0_deg),
            dsun_cm=float(dsun_cm),
            source="saved_observer_metadata",
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    carr_lon_deg = _as_float(_metadata_lookup(model_metadata, "crln_obs"))
    carr_lat_deg = _as_float(_metadata_lookup(model_metadata, "crlt_obs"))
    if carr_lon_deg is not None and carr_lat_deg is not None and dsun_cm is not None:
        converted = _carrington_observer_to_stonyhurst(carr_lon_deg, carr_lat_deg, dsun_cm, model_time)
        if converted is not None:
            return MetadataObserverState(
                observer_name=observer_name,
                l0_deg=converted[0],
                b0_deg=converted[1],
                dsun_cm=float(dsun_cm),
                source="model_metadata_carrington",
                rsun_cm=rsun_cm,
                rsun_arcsec=rsun_arcsec,
            )

    return None


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
    metadata_observer = _metadata_observer_state(model_metadata, observer_metadata, model_time)

    cli_name = normalize_observer_name(getattr(cli_args, "observer", None))
    if cli_name:
        try:
            l0_deg, b0_deg, dsun_cm = _observer_from_sunpy(cli_name, model_time)
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
        except Exception as exc:
            warnings.append(f"CLI observer '{getattr(cli_args, 'observer', cli_name)}' could not be resolved: {exc}")

    any_cli_triad = any(getattr(cli_args, key, None) is not None for key in ("lonc_deg", "b0sun_deg", "dsun_cm"))
    cli_lonc, cli_b0, cli_dsun = _cli_render_triad(cli_args, model)
    if any_cli_triad:
        actual_name = "earth"
        rsun_cm = None
        rsun_arcsec = None
        if metadata_observer is not None:
            l0_deg = metadata_observer.l0_deg
            b0_deg = metadata_observer.b0_deg
            dsun_cm = metadata_observer.dsun_cm
            actual_name = metadata_observer.observer_name
            rsun_cm = metadata_observer.rsun_cm
            rsun_arcsec = metadata_observer.rsun_arcsec
        else:
            meta_name = normalize_observer_name(
                _nested_lookup(observer_metadata, "name")
                or _metadata_lookup(model_metadata, "observer_name", "observer", "observatory", "obsrvtry")
            )
            if meta_name:
                try:
                    l0_deg, b0_deg, dsun_cm = _observer_from_sunpy(meta_name, model_time)
                    actual_name = meta_name
                except Exception as exc:
                    warnings.append(f"Model observer '{meta_name}' could not be resolved: {exc}")
                    l0_deg, b0_deg, dsun_cm = _observer_from_sunpy("earth", model_time)
            else:
                l0_deg, b0_deg, dsun_cm = _observer_from_sunpy("earth", model_time)
        rsun_cm, rsun_arcsec = _finalize_rsun_state(
            dsun_cm=dsun_cm,
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )
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

    if metadata_observer is not None:
        render_lonc = _render_lonc_for_observer(
            model_metadata,
            observer_lon_deg=metadata_observer.l0_deg,
            observer_lat_deg=metadata_observer.b0_deg,
            observer_dsun_cm=metadata_observer.dsun_cm,
            model_time=model_time,
            fallback_lonc_deg=model_lonc,
        )
        rsun_cm, rsun_arcsec = _finalize_rsun_state(
            dsun_cm=metadata_observer.dsun_cm,
            rsun_cm=metadata_observer.rsun_cm,
            rsun_arcsec=metadata_observer.rsun_arcsec,
        )
        return ResolvedObserverGeometry(
            metadata_observer.observer_name,
            metadata_observer.l0_deg,
            metadata_observer.b0_deg,
            metadata_observer.dsun_cm,
            render_lonc,
            metadata_observer.b0_deg,
            metadata_observer.dsun_cm,
            metadata_observer.source,
            tuple(warnings),
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )

    meta_name = normalize_observer_name(
        _nested_lookup(observer_metadata, "name")
        or _metadata_lookup(model_metadata, "observer_name", "observer", "observatory", "obsrvtry")
    )
    if meta_name:
        try:
            l0_deg, b0_deg, dsun_cm = _observer_from_sunpy(meta_name, model_time)
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
        except Exception as exc:
            warnings.append(f"Model observer '{meta_name}' could not be resolved: {exc}")

    l0_deg, b0_deg, dsun_cm = _observer_from_sunpy("earth", model_time)
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


def _metadata_square_fov(observer_metadata: dict[str, Any] | None) -> bool:
    for path in (("fov", "square"), ("fov_box", "square")):
        value = _nested_lookup(observer_metadata, *path)
        if value is not None:
            return bool(value)
    return False


def _index_box_corners_hgs(
    model: Any,
    model_metadata: dict[str, Any] | None,
    *,
    obstime: Time,
) -> SkyCoord | None:
    if not isinstance(model_metadata, dict):
        return None
    lon_ref = _as_float(_metadata_lookup(model_metadata, "crval1", "lon"))
    lat_ref = _as_float(_metadata_lookup(model_metadata, "crval2", "lat"))
    dsun_obs_m = _as_float(_metadata_lookup(model_metadata, "dsun_obs"))
    hgln_obs = _as_float(_metadata_lookup(model_metadata, "hgln_obs"))
    hglt_obs = _as_float(_metadata_lookup(model_metadata, "hglt_obs", "solar_b0"))
    rsun_ref_m = _as_float(_metadata_lookup(model_metadata, "rsun_ref"))
    box_nx = _as_float(_metadata_lookup(model_metadata, "box_nx"))
    box_ny = _as_float(_metadata_lookup(model_metadata, "box_ny"))
    box_nz = _as_float(_metadata_lookup(model_metadata, "box_nz"))
    box_dr_x = _as_float(_metadata_lookup(model_metadata, "box_dr_x"))
    box_dr_y = _as_float(_metadata_lookup(model_metadata, "box_dr_y"))
    box_dr_z = _as_float(_metadata_lookup(model_metadata, "box_dr_z"))
    if None in (lon_ref, lat_ref, dsun_obs_m, hgln_obs, hglt_obs, rsun_ref_m, box_nx, box_ny, box_nz, box_dr_x, box_dr_y, box_dr_z):
        return None

    try:
        anchor_hgs = SkyCoord(
            lon=float(lon_ref) * u.deg,
            lat=float(lat_ref) * u.deg,
            radius=float(rsun_ref_m) * u.m,
            frame=frames.HeliographicCarrington,
            observer="self",
            obstime=obstime,
        ).transform_to(frames.HeliographicStonyhurst(obstime=obstime))
    except Exception:
        return None

    center = np.array(
        [
            anchor_hgs.cartesian.x.to_value(u.cm),
            anchor_hgs.cartesian.y.to_value(u.cm),
            anchor_hgs.cartesian.z.to_value(u.cm),
        ],
        dtype=np.float64,
    )
    radial = center / float(np.linalg.norm(center))
    north_ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    y_axis = north_ref - np.dot(north_ref, radial) * radial
    y_norm = float(np.linalg.norm(y_axis))
    if not np.isfinite(y_norm) or y_norm <= 0:
        east_ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        y_axis = east_ref - np.dot(east_ref, radial) * radial
        y_norm = float(np.linalg.norm(y_axis))
        if not np.isfinite(y_norm) or y_norm <= 0:
            return None
    y_axis /= y_norm
    x_axis = np.cross(y_axis, radial)
    x_axis /= float(np.linalg.norm(x_axis))

    rsun_cm = float(model["RSun"][0])
    xdim_cm = float(box_nx) * float(box_dr_x) * rsun_cm
    ydim_cm = float(box_ny) * float(box_dr_y) * rsun_cm
    zdim_cm = float(box_nz) * float(box_dr_z) * rsun_cm

    rows = []
    for z_cm in (0.0, zdim_cm):
        for y_cm in (-0.5 * ydim_cm, 0.5 * ydim_cm):
            for x_cm in (-0.5 * xdim_cm, 0.5 * xdim_cm):
                rows.append(center + x_cm * x_axis + y_cm * y_axis + z_cm * radial)
    rows = np.asarray(rows, dtype=np.float64)
    return SkyCoord(
        x=rows[:, 0] * u.cm,
        y=rows[:, 1] * u.cm,
        z=rows[:, 2] * u.cm,
        representation_type="cartesian",
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
    )


def _execute_box_corners_hgs(
    model_metadata: dict[str, Any] | None,
    *,
    obstime: Time,
) -> SkyCoord | None:
    execute_text = _metadata_lookup(model_metadata, "execute")
    if not isinstance(execute_text, str):
        return None
    geometry = extract_geometry_from_execute(execute_text)
    if not isinstance(geometry, dict):
        return None

    observer_name = normalize_observer_name(geometry.get("geometry_observer")) or "earth"
    try:
        obs_lon_deg, obs_lat_deg, obs_dsun_cm = _observer_from_sunpy(observer_name, obstime)
    except Exception:
        if observer_name != "earth":
            obs_lon_deg, obs_lat_deg, obs_dsun_cm = _observer_from_sunpy("earth", obstime)
        else:
            return None

    geometry_observer = SkyCoord(
        lon=float(obs_lon_deg) * u.deg,
        lat=float(obs_lat_deg) * u.deg,
        radius=float(obs_dsun_cm) * u.cm,
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
    )
    rsun_cm = _RSUN_METERS * 100.0
    coord_mode = str(geometry.get("coord_mode") or "hpc").strip().lower()
    center_x = float(geometry["center_x"])
    center_y = float(geometry["center_y"])
    if coord_mode == "hpc":
        box_origin = SkyCoord(
            Tx=center_x * u.arcsec,
            Ty=center_y * u.arcsec,
            obstime=obstime,
            observer=geometry_observer,
            frame=frames.Helioprojective,
        )
    elif coord_mode == "hgc":
        box_origin = SkyCoord(
            lon=center_x * u.deg,
            lat=center_y * u.deg,
            radius=rsun_cm * u.cm,
            obstime=obstime,
            observer=geometry_observer,
            frame=frames.HeliographicCarrington,
        )
    else:
        box_origin = SkyCoord(
            lon=center_x * u.deg,
            lat=center_y * u.deg,
            radius=rsun_cm * u.cm,
            obstime=obstime,
            observer=geometry_observer,
            frame=frames.HeliographicStonyhurst,
        )

    box_origin_hgs = box_origin.transform_to(frames.HeliographicStonyhurst(obstime=obstime))
    local_frame = frames.Heliocentric(observer=box_origin_hgs, obstime=obstime)
    if hasattr(frames.Helioprojective, "assume_spherical_screen"):
        screen_ctx = frames.Helioprojective.assume_spherical_screen(geometry_observer)
    elif SphericalScreen is not None:
        screen_ctx = SphericalScreen(geometry_observer)
    else:  # pragma: no cover
        from contextlib import nullcontext
        screen_ctx = nullcontext()
    with screen_ctx:
        origin_local = box_origin_hgs.transform_to(local_frame)

    dims = tuple(int(v) for v in geometry["dims"])
    dim_mm = (np.asarray(dims, dtype=np.float64) * float(geometry["dx_km"])) / 1000.0
    center_local = SkyCoord(
        x=origin_local.x,
        y=origin_local.y,
        z=origin_local.z + (0.5 * dim_mm[2]) * u.Mm,
        frame=origin_local.frame,
    )
    half_x = 0.5 * dim_mm[0] * u.Mm
    half_y = 0.5 * dim_mm[1] * u.Mm
    half_z = 0.5 * dim_mm[2] * u.Mm
    corners_local = SkyCoord(
        x=[
            center_local.x - half_x, center_local.x + half_x, center_local.x - half_x, center_local.x + half_x,
            center_local.x - half_x, center_local.x + half_x, center_local.x - half_x, center_local.x + half_x,
        ],
        y=[
            center_local.y - half_y, center_local.y - half_y, center_local.y + half_y, center_local.y + half_y,
            center_local.y - half_y, center_local.y - half_y, center_local.y + half_y, center_local.y + half_y,
        ],
        z=[
            center_local.z - half_z, center_local.z - half_z, center_local.z - half_z, center_local.z - half_z,
            center_local.z + half_z, center_local.z + half_z, center_local.z + half_z, center_local.z + half_z,
        ],
        frame=center_local.frame,
    )
    return corners_local.transform_to(frames.HeliographicStonyhurst(obstime=obstime))


def _model_box_corners_hgs(model: Any) -> SkyCoord:
    obstime = model_time_from_model(model)
    lon = float(model["lonC"][0])
    lat = float(model["latC"][0])
    rsun_cm = float(model["RSun"][0])
    nx = int(model["Nx"][0])
    ny = int(model["Ny"][0])
    dx_cm = float(model["dx"][0])
    dy_cm = float(model["dy"][0])

    dz = np.asarray(model["dz"][0], dtype=np.float64)
    if dz.ndim != 3:
        raise ValueError(f"Unexpected model dz shape: {dz.shape}")
    z_top_cm = float(np.max(np.sum(dz, axis=0)))

    half_x_cm = 0.5 * float(nx) * dx_cm
    half_y_cm = 0.5 * float(ny) * dy_cm

    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    radial = np.array(
        [
            np.cos(lat_rad) * np.cos(lon_rad),
            np.cos(lat_rad) * np.sin(lon_rad),
            np.sin(lat_rad),
        ],
        dtype=np.float64,
    )
    center = radial * rsun_cm
    earth_obs = get_body_heliographic_stonyhurst("earth", obstime)
    earth_pos = np.array(
        [
            earth_obs.cartesian.x.to_value(u.cm),
            earth_obs.cartesian.y.to_value(u.cm),
            earth_obs.cartesian.z.to_value(u.cm),
        ],
        dtype=np.float64,
    )
    los = earth_pos - center
    los_norm = float(np.linalg.norm(los))
    if not np.isfinite(los_norm) or los_norm <= 0:
        raise ValueError("Invalid Earth LOS vector for model box construction.")
    los /= los_norm

    x_axis = np.cross(los, radial)
    x_norm = float(np.linalg.norm(x_axis))
    if not np.isfinite(x_norm) or x_norm <= 0:
        north_ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        y_axis = north_ref - np.dot(north_ref, radial) * radial
        y_norm = float(np.linalg.norm(y_axis))
        if not np.isfinite(y_norm) or y_norm <= 0:
            raise ValueError("Failed to construct tangent-plane basis.")
        y_axis /= y_norm
        x_axis = np.cross(y_axis, radial)
        x_axis /= float(np.linalg.norm(x_axis))
    else:
        x_axis /= x_norm
        y_axis = np.cross(radial, x_axis)
        y_axis /= float(np.linalg.norm(y_axis))

    cart_rows = []
    for z_cm in (0.0, z_top_cm):
        for y_cm in (-half_y_cm, half_y_cm):
            for x_cm in (-half_x_cm, half_x_cm):
                cart_rows.append(center + x_cm * x_axis + y_cm * y_axis + z_cm * radial)
    rows = np.asarray(cart_rows, dtype=np.float64)
    return SkyCoord(
        x=rows[:, 0] * u.cm,
        y=rows[:, 1] * u.cm,
        z=rows[:, 2] * u.cm,
        representation_type="cartesian",
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
    )


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
    hpc_frame = frames.Helioprojective(observer=observer, obstime=obstime)
    corners_hgs = _index_box_corners_hgs(model, model_metadata, obstime=obstime)
    if corners_hgs is None:
        corners_hgs = _execute_box_corners_hgs(model_metadata, obstime=obstime)
    if corners_hgs is None:
        corners_hgs = _model_box_corners_hgs(model)
    corners_hpc = corners_hgs.transform_to(hpc_frame)
    tx = np.asarray(corners_hpc.Tx.to_value(u.arcsec), dtype=np.float64)
    ty = np.asarray(corners_hpc.Ty.to_value(u.arcsec), dtype=np.float64)
    if tx.size != 8 or ty.size != 8 or not np.all(np.isfinite(tx)) or not np.all(np.isfinite(ty)):
        raise ValueError("Failed to project model box corners into observer frame.")
    pad = max(float(pad_arcsec), 0.0)
    xmin = float(np.min(tx)) - pad
    xmax = float(np.max(tx)) + pad
    ymin = float(np.min(ty)) - pad
    ymax = float(np.max(ty)) + pad
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
        "corners_hpc": corners_hpc,
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
    footprint = compute_inscribing_fov(
        model,
        observer_geometry,
        model_metadata=model_metadata,
        observer_metadata=observer_metadata,
        pad_arcsec=pad_xy_arcsec,
    )
    obstime = model_time_from_model(model)
    observer = build_observer_coordinate(observer_geometry, obstime)
    observer_frame = frames.Heliocentric(observer=observer, obstime=obstime)
    corners_hgs = _index_box_corners_hgs(model, model_metadata, obstime=obstime)
    if corners_hgs is None:
        corners_hgs = _execute_box_corners_hgs(model_metadata, obstime=obstime)
    if corners_hgs is None:
        corners_hgs = _model_box_corners_hgs(model)
    corners_hcc = corners_hgs.transform_to(observer_frame)
    z_vals = np.asarray(corners_hcc.z.to_value(u.Mm), dtype=np.float64)
    finite = np.isfinite(z_vals)
    if not np.any(finite):
        raise ValueError("Failed to derive observer-heliocentric z range.")
    z_vals = z_vals[finite]
    z_min = float(np.nanmin(z_vals))
    z_max = float(np.nanmax(z_vals))
    z_span = max(1e-6, z_max - z_min)
    z_pad = max(0.0, float(pad_z_frac)) * z_span
    result = dict(footprint)
    result.update(
        {
            "zmin_mm": z_min - z_pad,
            "zmax_mm": z_max + z_pad,
        }
    )
    return result


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
    ref_coord = SkyCoord(
        Tx=float(xc_arcsec) * u.arcsec,
        Ty=float(yc_arcsec) * u.arcsec,
        frame=frames.Helioprojective(observer=observer, obstime=observer.obstime),
    )
    header = make_fitswcs_header(
        np.empty((int(ny), int(nx)), dtype=np.float32),
        ref_coord,
        scale=u.Quantity([float(dx_arcsec), float(dy_arcsec)], u.arcsec / u.pix),
    )
    rsun_ref_m = (
        float(observer_geometry.rsun_cm) / 100.0
        if observer_geometry.rsun_cm is not None and np.isfinite(observer_geometry.rsun_cm)
        else _RSUN_METERS
    )
    rsun_obs_arcsec = observer_geometry.rsun_arcsec
    if rsun_obs_arcsec is None or not np.isfinite(rsun_obs_arcsec) or rsun_obs_arcsec <= 0:
        ratio = min(1.0, rsun_ref_m / (float(observer_geometry.dsun_cm) / 100.0))
        rsun_obs_arcsec = float(np.degrees(np.arcsin(ratio)) * 3600.0)
    header["DATE-OBS"] = Time(obs_time).isot
    header["BUNIT"] = str(bunit)
    header["OBSERVER"] = _pretty_observer_name(observer_geometry.observer_name)
    # Keep explicit observer-triad aliases alongside the standard WCS cards so
    # downstream readers can recover the saved LOS directly from each map.
    header["B0"] = float(observer_geometry.b0_deg)
    header["L0"] = float(observer_geometry.l0_deg)
    header["RSUN_ARC"] = float(rsun_obs_arcsec)
    header["SOLAR_B0"] = float(observer_geometry.b0_deg)
    header["SOLAR_L0"] = float(observer_geometry.l0_deg)
    header["HGLN_OBS"] = float(observer_geometry.l0_deg)
    header["HGLT_OBS"] = float(observer_geometry.b0_deg)
    header["DSUN_OBS"] = float(observer_geometry.dsun_cm) / 100.0
    header["RSUN_REF"] = float(rsun_ref_m)
    header["RSUN_OBS"] = float(rsun_obs_arcsec)
    try:
        observer_hgc = observer.transform_to(frames.HeliographicCarrington(obstime=observer.obstime, observer="self"))
        header["CRLN_OBS"] = float(observer_hgc.lon.to_value(u.deg))
        header["CRLT_OBS"] = float(observer_hgc.lat.to_value(u.deg))
    except Exception:
        pass
    return header
