from __future__ import annotations

import re
import shlex
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict

import astropy.units as u
import h5py
import numpy as np
import scipy.io as io
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from pyampp.geometry import project_world_to_observer_hpc
from pyampp.io import load_model
from sunpy.coordinates import frames, get_earth, sun

from .voxel_id import gx_box2id


_STANDARD_OBSERVER_ALIASES = {
    "earth": "earth",
    "terra": "earth",
    "sdo": "earth",
    "solar dynamics observatory": "earth",
    "solo": "solar orbiter",
    "solar orbiter": "solar orbiter",
    "solar-orbiter": "solar orbiter",
    "stereo-a": "stereo-a",
    "stereo a": "stereo-a",
    "stereoa": "stereo-a",
    "stereo-b": "stereo-b",
    "stereo b": "stereo-b",
    "stereob": "stereo-b",
}


def decode_if_bytes(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    return value


def _tokenize_execute(execute_text: str) -> list[str]:
    try:
        return [str(tok) for tok in shlex.split(str(execute_text))]
    except Exception:
        return [str(tok) for tok in str(execute_text).split()]


def _extract_from_execute(execute_text: str):
    tokens = _tokenize_execute(execute_text)
    dims = None
    dx_km = None
    box_res_km = None

    for idx, token in enumerate(tokens):
        lowered = token.strip().lower()
        if lowered == "--box-dims" and idx + 3 < len(tokens):
            try:
                dims = (int(tokens[idx + 1]), int(tokens[idx + 2]), int(tokens[idx + 3]))
            except ValueError:
                pass
        elif lowered.startswith("--box-dims="):
            parts = lowered.split("=", 1)[1].replace(",", " ").split()
            if len(parts) == 3:
                try:
                    dims = tuple(int(part) for part in parts)
                except ValueError:
                    pass
        elif lowered == "--dx-km" and idx + 1 < len(tokens):
            try:
                dx_km = float(tokens[idx + 1])
            except ValueError:
                pass
        elif lowered.startswith("--dx-km="):
            try:
                dx_km = float(token.split("=", 1)[1])
            except ValueError:
                pass
        elif lowered == "--box-res" and idx + 1 < len(tokens):
            try:
                box_res_km = float(tokens[idx + 1]) * 1000.0
            except ValueError:
                pass
        elif lowered.startswith("--box-res="):
            try:
                box_res_km = float(token.split("=", 1)[1]) * 1000.0
            except ValueError:
                pass

    dims_match = None
    dx_match = None
    box_res_match = None
    if dims_match is None:
        dims_match = re.search(
            r"SIZE_PIX\s*=\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]",
            execute_text,
            flags=re.IGNORECASE,
        )
    if dx_match is None:
        dx_match = re.search(r"DX_KM\s*=\s*([0-9.]+)", execute_text, flags=re.IGNORECASE)
    if box_res_match is None:
        box_res_match = re.search(r"BOX_RES\s*=\s*([0-9.]+)", execute_text, flags=re.IGNORECASE)
    if dims is None and dims_match:
        dims = (int(dims_match.group(1)), int(dims_match.group(2)), int(dims_match.group(3)))
    if dx_km is None and dx_match:
        dx_km = float(dx_match.group(1))
    elif dx_km is None and box_res_km is not None:
        dx_km = box_res_km
    elif dx_km is None and box_res_match:
        dx_km = float(box_res_match.group(1)) * 1000.0
    return dims, dx_km


def extract_geometry_from_execute(execute_text: str) -> dict[str, Any] | None:
    if not execute_text:
        return None
    tokens = _tokenize_execute(execute_text)

    center = extract_center_from_execute(str(execute_text))
    dims, dx_km = _extract_from_execute(str(execute_text))

    coord_mode = None
    geometry_observer = None
    projection = None
    for idx, token in enumerate(tokens):
        lowered = str(token).strip().lower()
        if lowered in {"--hpc", "--hgc", "--hgs"}:
            coord_mode = lowered[2:]
        elif lowered in {"--cea", "--tan"}:
            projection = lowered[2:]
        elif lowered == "--observer-name" and idx + 1 < len(tokens):
            geometry_observer = str(tokens[idx + 1]).strip()

    if coord_mode is None:
        if re.search(r"(?:^|\s)--hpc(?:\s|$)", execute_text):
            coord_mode = "hpc"
        elif re.search(r"(?:^|\s)--hgc(?:\s|$)", execute_text):
            coord_mode = "hgc"
        elif re.search(r"(?:^|\s)--hgs(?:\s|$)", execute_text):
            coord_mode = "hgs"

    if center is None or dims is None or dx_km is None:
        return None

    return {
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "coord_mode": str(coord_mode or "hpc"),
        "geometry_observer": str(geometry_observer or "earth"),
        "dims": tuple(int(v) for v in dims),
        "dx_km": float(dx_km),
        "projection": projection,
    }


def extract_center_from_execute(execute_text: str) -> tuple[float, float] | None:
    tokens = _tokenize_execute(execute_text)
    for idx, token in enumerate(tokens):
        lowered = token.strip().lower()
        if lowered == "--coords" and idx + 2 < len(tokens):
            try:
                return float(tokens[idx + 1]), float(tokens[idx + 2])
            except ValueError:
                pass
        if lowered.startswith("--coords="):
            parts = token.split("=", 1)[1].replace(",", " ").split()
            if len(parts) == 2:
                try:
                    return float(parts[0]), float(parts[1])
                except ValueError:
                    pass

    # Legacy IDL execute style: CENTER_ARCSEC=[ X, Y ]
    m_center = re.search(
        r"CENTER_ARCSEC\s*=\s*\[\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*\]",
        execute_text,
        flags=re.IGNORECASE,
    )
    if m_center:
        return float(m_center.group(1)), float(m_center.group(2))
    return None


def infer_center_from_execute(loader_name: str, model_path: Path) -> tuple[float, float] | None:
    execute_text = _load_execute_text_from_model(loader_name, model_path)

    if not execute_text:
        return None
    return extract_center_from_execute(str(execute_text))


def infer_fov_from_execute(
    loader_name: str,
    model_path: Path,
    dsun_cm: float,
    fallback_nx: int,
    fallback_ny: int,
    fallback_dx_cm: float,
) -> tuple[float, float]:
    execute_text = _load_execute_text_from_model(loader_name, model_path)

    dims, dx_km = (None, None)
    if execute_text:
        dims, dx_km = _extract_from_execute(execute_text)

    if dims is not None and dx_km is not None:
        model_w_km = dims[0] * dx_km
        model_h_km = dims[1] * dx_km
    else:
        model_w_km = fallback_nx * (fallback_dx_cm / 1e5)
        model_h_km = fallback_ny * (fallback_dx_cm / 1e5)

    dsun_km = dsun_cm / 1e5
    km_to_arcsec = 206265.0 / dsun_km
    fov_x = model_w_km * km_to_arcsec
    fov_y = model_h_km * km_to_arcsec
    return float(fov_x), float(fov_y)


def estimate_hpc_center(
    model,
    *,
    observer_lon_deg: float | None = None,
    observer_lat_deg: float | None = None,
    observer_dsun_cm: float | None = None,
) -> tuple[float, float]:
    unix = float(model["obstime"][0]) + 283996800.0
    obs_time = Time(unix, format="unix")
    if observer_lon_deg is None and observer_lat_deg is None and observer_dsun_cm is None:
        observer = get_earth(obs_time)
    else:
        dsun_cm = float(observer_dsun_cm if observer_dsun_cm is not None else sun.earth_distance(obs_time).to_value(u.cm))
        observer = SkyCoord(
            lon=float(observer_lon_deg if observer_lon_deg is not None else 0.0) * u.deg,
            lat=float(observer_lat_deg if observer_lat_deg is not None else 0.0) * u.deg,
            radius=dsun_cm * u.cm,
            frame=frames.HeliographicStonyhurst,
            obstime=obs_time,
        )
    center_hgs = SkyCoord(
        lon=float(model["lonC"][0]) * u.deg,
        lat=float(model["latC"][0]) * u.deg,
        radius=(float(model["RSun"][0]) / 1e5) * u.km,
        frame=frames.HeliographicStonyhurst,
        obstime=obs_time,
        observer=observer,
    )
    center_hpc = project_world_to_observer_hpc(center_hgs, observer=observer, obstime=obs_time)
    if center_hpc is None:
        center_hpc = center_hgs.transform_to(frames.Helioprojective(observer=observer, obstime=obs_time))
    return float(center_hpc.Tx.to_value(u.arcsec)), float(center_hpc.Ty.to_value(u.arcsec))


def _load_execute_text_from_model(loader_name: str, model_path: Path) -> str | None:
    try:
        model_payload = load_model(model_path, strict=False)
    except Exception:
        return None

    metadata = model_payload.get("metadata") if isinstance(model_payload, dict) else None
    if not isinstance(metadata, dict):
        return None
    execute_text = decode_if_bytes(metadata.get("execute"))
    if execute_text is None:
        return None
    return str(execute_text)


@dataclass
class ChromoModelData:
    header: Dict[str, Any]
    model: Dict[str, np.ndarray]
    observer: Any | None = None


def _coerce_time(obs_time):
    if isinstance(obs_time, Time):
        return obs_time
    if isinstance(obs_time, (bytes, np.bytes_)):
        obs_time = obs_time.decode("utf-8")
    return Time(obs_time)


def _normalize_observer_name(name: Any) -> str | None:
    value = decode_if_bytes(_scalar_from_any(name))
    if value is None:
        return None
    cleaned = " ".join(str(value).strip().split())
    if not cleaned:
        return None
    return _STANDARD_OBSERVER_ALIASES.get(cleaned.lower(), cleaned.lower())


def _file_observer_state(header: Dict[str, Any]) -> dict[str, Any]:
    state: dict[str, Any] = {}
    observer_name = _normalize_observer_name(header.get("observer_name", header.get("observer")))
    if observer_name:
        state["observer_name"] = observer_name

    dsun_cm = header.get("observer_dsun_cm", None)
    if dsun_cm is None:
        dsun_cm = header.get("dsun_obs", None)
    if dsun_cm is not None:
        dsun_cm = float(dsun_cm)
        if dsun_cm < 1e12:
            dsun_cm *= 100.0
        state["DSun"] = dsun_cm

    b0sun = header.get("observer_b0_deg", None)
    if b0sun is None:
        b0sun = header.get("solar_b0", None)
    if b0sun is None:
        b0sun = header.get("crlt_obs", None)
    if b0sun is None:
        b0sun = header.get("hglt_obs", None)
    if b0sun is not None:
        state["b0Sun"] = float(b0sun)

    lon_deg = header.get("observer_hgln_obs_deg", None)
    if lon_deg is None:
        lon_deg = header.get("hgln_obs", None)
    if lon_deg is not None:
        state["observer_lon_deg"] = float(lon_deg)

    crlon_deg = header.get("crln_obs", None)
    if crlon_deg is None:
        crlon_deg = header.get("observer_l0_deg", None)
    if crlon_deg is not None:
        state["observer_crln_obs_deg"] = float(crlon_deg)

    lat_deg = header.get("observer_hglt_obs_deg", None)
    if lat_deg is None:
        lat_deg = header.get("hglt_obs", None)
    if lat_deg is None:
        lat_deg = header.get("crlt_obs", None)
    if lat_deg is not None:
        state["observer_lat_deg"] = float(lat_deg)

    return state


def _recompute_observer_state(header: Dict[str, Any], obs_time: Time, observer_name: str | None) -> dict[str, Any]:
    state: dict[str, Any] = {}
    normalized = _normalize_observer_name(observer_name) or "earth"
    state["observer_name"] = normalized

    if normalized == "earth":
        observer = get_earth(obs_time)
        state["observer_crln_obs_deg"] = float(sun.L0(obs_time).to(u.deg).value)
    elif normalized == "solar orbiter":
        from sunpy.coordinates.ephemeris import get_horizons_coord

        observer = get_horizons_coord("Solar Orbiter", obs_time)
    elif normalized == "stereo-a":
        from sunpy.coordinates.ephemeris import get_horizons_coord

        observer = get_horizons_coord("STEREO-A", obs_time)
    elif normalized == "stereo-b":
        from sunpy.coordinates.ephemeris import get_horizons_coord

        observer = get_horizons_coord("STEREO-B", obs_time)
    else:
        raise ValueError(f"Cannot recompute observer ephemeris for observer '{observer_name}'.")

    state["DSun"] = float(observer.radius.to(u.cm).value)
    state["b0Sun"] = float(observer.lat.to(u.deg).value)
    state["observer_lon_deg"] = float(observer.lon.to(u.deg).value)
    state["observer_lat_deg"] = float(observer.lat.to(u.deg).value)
    if "observer_crln_obs_deg" not in state:
        try:
            carr = observer.transform_to(frames.HeliographicCarrington(observer="self", obstime=obs_time))
            state["observer_crln_obs_deg"] = float(carr.lon.to(u.deg).value)
        except Exception:
            pass
    return state


def _sanitize_status_mask(start_idx, end_idx, chromo_mask):
    max_index = int(np.prod(chromo_mask.shape))
    start_idx = np.asarray(start_idx, dtype=np.int64).copy()
    end_idx = np.asarray(end_idx, dtype=np.int64).copy()
    start_idx = np.clip(start_idx, 0, max_index - 1)
    end_idx = np.clip(end_idx, 0, max_index - 1)
    return start_idx, end_idx


def _load_h5_group(group: h5py.Group) -> Dict[str, np.ndarray]:
    loaded = {}
    for key in group.keys():
        try:
            loaded[key] = group[key][:]
        except ValueError:
            loaded[key] = group[key][()]
    return loaded


def _normalize_cube_xyz(arr: np.ndarray, nx_hint: int | None = None, ny_hint: int | None = None) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim != 3:
        return a

    if nx_hint is not None and ny_hint is not None:
        if a.shape[0] == nx_hint and a.shape[1] == ny_hint:
            return a
        if a.shape[0] == ny_hint and a.shape[1] == nx_hint:
            return a.transpose((1, 0, 2))
        if a.shape[1] == ny_hint and a.shape[2] == nx_hint:
            return a.transpose((2, 1, 0))
        if a.shape[1] == nx_hint and a.shape[2] == ny_hint:
            return a.transpose((1, 2, 0))

    return a.transpose((1, 0, 2))


def _normalize_vector_cube_xyzc(arr: np.ndarray, nx_hint: int | None = None, ny_hint: int | None = None) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim != 4:
        return a

    if a.shape[-1] == 3:
        comps = [_normalize_cube_xyz(a[..., i], nx_hint, ny_hint) for i in range(3)]
        return np.stack(comps, axis=-1)
    if a.shape[0] == 3:
        comps = [_normalize_cube_xyz(a[i, ...], nx_hint, ny_hint) for i in range(3)]
        return np.stack(comps, axis=-1)
    return a


def _components_to_vector_cube_xyzc(
    bx: np.ndarray, by: np.ndarray, bz: np.ndarray, nx_hint: int | None = None, ny_hint: int | None = None
) -> np.ndarray:
    bx_i = _normalize_cube_xyz(bx, nx_hint, ny_hint)
    by_i = _normalize_cube_xyz(by, nx_hint, ny_hint)
    bz_i = _normalize_cube_xyz(bz, nx_hint, ny_hint)
    if bx_i.shape != by_i.shape or bx_i.shape != bz_i.shape:
        raise ValueError(f"Incompatible vector component shapes: bx={bx_i.shape}, by={by_i.shape}, bz={bz_i.shape}")
    return np.stack((bx_i, by_i, bz_i), axis=-1)


def _scalar_from_any(value: Any) -> Any:
    arr = np.asarray(value)
    if arr.ndim == 0:
        return arr.item()
    if arr.size == 1:
        return arr.reshape(-1)[0].item()
    return value


def _flatten_status_cube(arr: np.ndarray) -> np.ndarray:
    a = np.asarray(arr)
    if a.ndim <= 1:
        return a.reshape(-1)
    return a.T.reshape(-1, order="F")


def _parse_execute_time(execute_text: str) -> str | None:
    try:
        tokens = shlex.split(execute_text)
    except ValueError:
        tokens = execute_text.split()
    for i, tok in enumerate(tokens):
        if tok == "--time" and i + 1 < len(tokens):
            return tokens[i + 1]
        if tok.startswith("--time="):
            return tok.split("=", 1)[1]
    return None


def _decode_dataset_scalar(ds) -> Any:
    value = ds[()]
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", "ignore")
    return value


def _extract_observer_group_metadata(model_f: h5py.File, header: Dict[str, Any]) -> None:
    observer_group = model_f.get("observer")
    if not isinstance(observer_group, h5py.Group):
        return

    def _copy_scalar(group: h5py.Group, source_key: str, target_key: str) -> None:
        if target_key in header or source_key not in group:
            return
        try:
            header[target_key] = decode_if_bytes(_decode_dataset_scalar(group[source_key]))
        except Exception:
            return

    _copy_scalar(observer_group, "name", "observer_name")
    _copy_scalar(observer_group, "label", "observer_label")
    _copy_scalar(observer_group, "source", "observer_source")

    ephemeris = observer_group.get("ephemeris")
    if isinstance(ephemeris, h5py.Group):
        for source_key, target_key in (
            ("obs_date", "observer_obs_date"),
            ("obs_time", "observer_obs_time"),
            ("hgln_obs_deg", "observer_hgln_obs_deg"),
            ("hglt_obs_deg", "observer_hglt_obs_deg"),
            ("dsun_cm", "observer_dsun_cm"),
            ("rsun_cm", "observer_rsun_cm"),
        ):
            _copy_scalar(ephemeris, source_key, target_key)

    pb0r = observer_group.get("pb0r")
    if isinstance(pb0r, h5py.Group):
        for source_key, target_key in (
            ("obs_date", "observer_pb0r_obs_date"),
            ("b0_deg", "observer_b0_deg"),
            ("l0_deg", "observer_l0_deg"),
            ("p_deg", "observer_p_deg"),
            ("rsun_arcsec", "observer_rsun_arcsec"),
        ):
            _copy_scalar(pb0r, source_key, target_key)

    fov = observer_group.get("fov")
    if isinstance(fov, h5py.Group):
        for source_key, target_key in (
            ("xc_arcsec", "observer_fov_xc_arcsec"),
            ("yc_arcsec", "observer_fov_yc_arcsec"),
            ("xsize_arcsec", "observer_fov_xsize_arcsec"),
            ("ysize_arcsec", "observer_fov_ysize_arcsec"),
            ("frame", "observer_fov_frame"),
            ("square", "observer_fov_square"),
        ):
            _copy_scalar(fov, source_key, target_key)

    if "observer" not in header and "observer_name" in header:
        header["observer"] = header["observer_name"]


def _decode_h5_group_raw(group: h5py.Group) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in group.items():
        if isinstance(value, h5py.Group):
            result[key] = _decode_h5_group_raw(value)
        else:
            raw = value[()]
            if isinstance(raw, (bytes, np.bytes_)):
                result[key] = raw.decode("utf-8", "ignore")
            elif isinstance(raw, np.ndarray) and raw.shape == ():
                item = raw.item()
                if isinstance(item, (bytes, np.bytes_)):
                    result[key] = item.decode("utf-8", "ignore")
                else:
                    result[key] = item
            else:
                result[key] = raw
    return result


def _base_index_fits_header_from_value(raw_index: Any) -> fits.Header | None:
    if raw_index is None:
        return None
    text = decode_if_bytes(raw_index)
    normalized = str(text).replace("\r\n", "\n").replace("\r", "\n").strip()
    if not normalized:
        return None
    first_line = normalized.split("\n", 1)[0].lstrip()
    # Older parity clone fixtures stored ``str(index)`` for the full structured
    # SAV record here rather than FITS header text. Skip those tuple repr blobs
    # instead of passing them to Astropy's FITS card parser.
    if first_line.startswith("("):
        return None
    try:
        header = fits.Header.fromstring(normalized, sep="\n")
    except Exception:
        return None
    return header if len(header) > 0 else None


def _header_float_value(header: fits.Header, *keys: str) -> float | None:
    for key in keys:
        value = header.get(key)
        if value is None:
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def _fill_header_from_base_index(base_index_raw: Any, header: Dict[str, Any]) -> None:
    # Handle both raw index data and HDF5 file objects
    if hasattr(base_index_raw, "get") and "base" in base_index_raw:
        # This is an HDF5 file object with a "base" group
        base_box = base_index_raw.get("base")
        if isinstance(base_box, dict):
            base_index_raw = base_box.get("index")
        else:
            # Handle h5py Group objects
            base_index_raw = base_index_raw["base"]["index"][()] if "base" in base_index_raw and "index" in base_index_raw["base"] else None

    base_index = _base_index_fits_header_from_value(base_index_raw)
    if base_index is None:
        return

    if "lon" not in header:
        lon = _header_float_value(base_index, "CRVAL1")
        if lon is not None:
            header["lon"] = lon
    if "lat" not in header:
        lat = _header_float_value(base_index, "CRVAL2")
        if lat is not None:
            header["lat"] = lat
    if "obs_time" not in header:
        date_obs = base_index.get("DATE-OBS", base_index.get("DATE_OBS"))
        if date_obs:
            header["obs_time"] = str(date_obs)
    if "dsun_obs" not in header:
        dsun_obs = _header_float_value(base_index, "DSUN_OBS")
        if dsun_obs is not None:
            header["dsun_obs"] = dsun_obs
    if "observer" not in header:
        observer_name = base_index.get("OBSERVER", base_index.get("OBSERVATORY"))
        if observer_name:
            header["observer"] = str(observer_name)
    if "hgln_obs" not in header:
        hglon = _header_float_value(base_index, "HGLN_OBS")
        if hglon is not None:
            header["hgln_obs"] = hglon
    if "hglt_obs" not in header:
        hglat = _header_float_value(base_index, "HGLT_OBS")
        if hglat is not None:
            header["hglt_obs"] = hglat
    if "crln_obs" not in header:
        crlon = _header_float_value(base_index, "CRLN_OBS")
        if crlon is not None:
            header["crln_obs"] = crlon
    if "crlt_obs" not in header:
        crlat = _header_float_value(base_index, "CRLT_OBS")
        if crlat is not None:
            header["crlt_obs"] = crlat
    # Do not infer lonC from index fallback tags (HGLN_OBS/CRVAL1).
    # Keep lonC derived from coordinate transforms in load_model_dict,
    # unless an explicit lonC override is provided by higher-level callers.


def _extract_observer_payload_metadata(observer_payload: Any, header: Dict[str, Any]) -> None:
    if not isinstance(observer_payload, dict):
        return

    def _copy_scalar(source: dict[str, Any], source_key: str, target_key: str) -> None:
        if target_key in header or source_key not in source:
            return
        header[target_key] = decode_if_bytes(_scalar_from_any(source[source_key]))

    _copy_scalar(observer_payload, "name", "observer_name")
    _copy_scalar(observer_payload, "label", "observer_label")
    _copy_scalar(observer_payload, "source", "observer_source")

    ephemeris = observer_payload.get("ephemeris")
    if isinstance(ephemeris, dict):
        for source_key, target_key in (
            ("obs_date", "observer_obs_date"),
            ("obs_time", "observer_obs_time"),
            ("hgln_obs_deg", "observer_hgln_obs_deg"),
            ("hglt_obs_deg", "observer_hglt_obs_deg"),
            ("dsun_cm", "observer_dsun_cm"),
            ("rsun_cm", "observer_rsun_cm"),
        ):
            _copy_scalar(ephemeris, source_key, target_key)

    pb0r = observer_payload.get("pb0r")
    if isinstance(pb0r, dict):
        for source_key, target_key in (
            ("obs_date", "observer_pb0r_obs_date"),
            ("b0_deg", "observer_b0_deg"),
            ("l0_deg", "observer_l0_deg"),
            ("p_deg", "observer_p_deg"),
            ("rsun_arcsec", "observer_rsun_arcsec"),
        ):
            _copy_scalar(pb0r, source_key, target_key)

    fov = observer_payload.get("fov")
    if isinstance(fov, dict):
        for source_key, target_key in (
            ("xc_arcsec", "observer_fov_xc_arcsec"),
            ("yc_arcsec", "observer_fov_yc_arcsec"),
            ("xsize_arcsec", "observer_fov_xsize_arcsec"),
            ("ysize_arcsec", "observer_fov_ysize_arcsec"),
            ("frame", "observer_fov_frame"),
            ("square", "observer_fov_square"),
        ):
            _copy_scalar(fov, source_key, target_key)

    if "observer" not in header and "observer_name" in header:
        header["observer"] = header["observer_name"]


def _geometry_contract_to_dict(contract: Any) -> dict[str, Any] | None:
    if contract is None:
        return None
    if isinstance(contract, dict):
        return contract
    to_dict = getattr(contract, "to_dict", None)
    if callable(to_dict):
        try:
            value = to_dict()
            if isinstance(value, dict):
                return value
        except Exception:
            return None
    return None


def _apply_geometry_contract_to_header(contract_dict: dict[str, Any], header: Dict[str, Any]) -> None:
    for src_key, dst_key in (
        ("anchor_lon_deg", "lon"),
        ("anchor_lat_deg", "lat"),
        ("obstime", "obs_time"),
        ("nx", "box_nx"),
        ("ny", "box_ny"),
        ("nz", "box_nz"),
        ("dr_x", "box_dr_x"),
        ("dr_y", "box_dr_y"),
        ("dr_z", "box_dr_z"),
    ):
        if dst_key in header or src_key not in contract_dict:
            continue
        header[dst_key] = decode_if_bytes(contract_dict[src_key])


def _normalize_observer_payload(observer_data: dict[str, Any] | None) -> dict[str, Any] | None:
    """Decode bytes values in observer payload for consistent string handling."""
    if not isinstance(observer_data, dict):
        return observer_data

    normalized = {}
    for key, value in observer_data.items():
        if isinstance(value, dict):
            # Recursively normalize nested dicts
            normalized[key] = _normalize_observer_payload(value)
        elif isinstance(value, (bytes, np.bytes_)):
            # Decode bytes to strings
            try:
                normalized[key] = value.decode('utf-8')
            except (AttributeError, UnicodeDecodeError):
                normalized[key] = value
        else:
            normalized[key] = value
    return normalized


def _build_chromo_data_from_payload(model_payload: Dict[str, Any]) -> ChromoModelData:
    chromo_box = model_payload.get("chromo")
    if not isinstance(chromo_box, dict):
        raise KeyError("Missing required '/chromo' payload in canonical model dictionary.")

    corona_box = model_payload.get("corona") if isinstance(model_payload.get("corona"), dict) else {}
    lines_box = model_payload.get("lines") if isinstance(model_payload.get("lines"), dict) else {}
    base_box = model_payload.get("base") if isinstance(model_payload.get("base"), dict) else {}
    metadata_box = model_payload.get("metadata") if isinstance(model_payload.get("metadata"), dict) else {}
    observer_raw = model_payload.get("observer") if isinstance(model_payload.get("observer"), dict) else None

    model_dict = {
        k: v for k, v in chromo_box.items() if k not in {"attrs", "bx", "by", "bz"}
    }
    header: Dict[str, Any] = {}

    chromo_attrs = chromo_box.get("attrs")
    if isinstance(chromo_attrs, dict):
        for key, value in chromo_attrs.items():
            header[key] = decode_if_bytes(value)

    corona_attrs = corona_box.get("attrs")
    if isinstance(corona_attrs, dict):
        for key, value in corona_attrs.items():
            if key not in header:
                header[key] = decode_if_bytes(value)

    for key, value in metadata_box.items():
        if key in header or key == "geometry_contract":
            continue
        arr = np.asarray(value)
        if arr.ndim > 0 and arr.size > 1:
            continue
        header[key] = decode_if_bytes(_scalar_from_any(value))

    contract_dict = _geometry_contract_to_dict(metadata_box.get("geometry_contract"))
    if "geometry_contract" in metadata_box and metadata_box.get("geometry_contract") is not None:
        header["geometry_contract"] = metadata_box.get("geometry_contract")
    if isinstance(contract_dict, dict):
        _apply_geometry_contract_to_header(contract_dict, header)

    ny_hint = None
    nx_hint = None
    base_bx = base_box.get("bx")
    if base_bx is not None:
        base_shape = np.asarray(base_bx).shape
        if len(base_shape) == 2:
            ny_hint, nx_hint = int(base_shape[0]), int(base_shape[1])

    axis_order_3d = None
    if "axis_order_3d" in metadata_box:
        try:
            axis_order_3d = str(decode_if_bytes(_scalar_from_any(metadata_box["axis_order_3d"]))).strip().lower()
        except Exception:
            axis_order_3d = None

    if "bx" in corona_box:
        raw_shape = np.asarray(corona_box["bx"]).shape
        if len(raw_shape) == 3:
            if axis_order_3d == "zyx":
                header.setdefault("box_nx", int(raw_shape[2]))
                header.setdefault("box_ny", int(raw_shape[1]))
                header.setdefault("box_nz", int(raw_shape[0]))
            else:
                header.setdefault("box_nx", int(raw_shape[0]))
                header.setdefault("box_ny", int(raw_shape[1]))
                header.setdefault("box_nz", int(raw_shape[2]))

    if "dr" not in model_dict and "dr" in corona_box:
        model_dict["dr"] = corona_box["dr"]
    if "dr" in corona_box:
        try:
            dr_arr = np.asarray(corona_box["dr"], dtype=float).reshape(-1)
            if dr_arr.size >= 3:
                header.setdefault("box_dr_x", float(dr_arr[0]))
                header.setdefault("box_dr_y", float(dr_arr[1]))
                header.setdefault("box_dr_z", float(dr_arr[2]))
        except Exception:
            pass
    if "dr" not in model_dict and "dr" in chromo_box:
        model_dict["dr"] = chromo_box["dr"]

    if "corona_base" not in model_dict:
        if "corona_base" in corona_box:
            model_dict["corona_base"] = corona_box["corona_base"]
        elif "corona_base" in chromo_box:
            model_dict["corona_base"] = chromo_box["corona_base"]

    if "bcube" not in model_dict:
        if all(k in corona_box for k in ("bx", "by", "bz")):
            model_dict["bcube"] = _components_to_vector_cube_xyzc(
                corona_box["bx"], corona_box["by"], corona_box["bz"], nx_hint, ny_hint
            )
        elif "bcube" in corona_box:
            model_dict["bcube"] = corona_box["bcube"]
        elif "bcube" in chromo_box:
            model_dict["bcube"] = chromo_box["bcube"]

    if "chromo_bcube" not in model_dict:
        if all(k in chromo_box for k in ("bx", "by", "bz")):
            model_dict["chromo_bcube"] = _components_to_vector_cube_xyzc(
                chromo_box["bx"], chromo_box["by"], chromo_box["bz"], nx_hint, ny_hint
            )
        elif "chromo_bcube" in chromo_box:
            model_dict["chromo_bcube"] = chromo_box["chromo_bcube"]

    if "chromo_mask" not in model_dict:
        if "chromo_mask" in chromo_box:
            model_dict["chromo_mask"] = chromo_box["chromo_mask"]
        elif "chromo_mask" in base_box:
            model_dict["chromo_mask"] = base_box["chromo_mask"]

    for key in ("av_field", "phys_length", "voxel_status", "start_idx", "end_idx"):
        if key not in model_dict and key in lines_box:
            model_dict[key] = lines_box[key]
        if key not in model_dict and key in chromo_box:
            model_dict[key] = chromo_box[key]

    _extract_observer_payload_metadata(observer_raw, header)
    _fill_header_from_base_index(base_box.get("index"), header)
    if "obs_time" not in header:
        for key in ("observer_pb0r_obs_date", "observer_obs_date", "observer_obs_time"):
            value = header.get(key)
            if value:
                header["obs_time"] = value
                break

    execute_text = header.get("execute")
    if "dr" not in model_dict and isinstance(execute_text, str):
        _dims, dx_km = _extract_from_execute(execute_text)
        if dx_km is not None:
            rsun_cm = float(header.get("observer_rsun_cm", 69600000000.0))
            rsun_km = rsun_cm / 1e5
            dr = float(dx_km) / rsun_km
            model_dict["dr"] = np.array([dr, dr, dr], dtype=np.float64)

    if "dz" in model_dict and np.asarray(model_dict["dz"]).ndim == 3:
        model_dict["dz"] = _normalize_cube_xyz(np.asarray(model_dict["dz"]), nx_hint, ny_hint)
    if "bcube" in model_dict and np.asarray(model_dict["bcube"]).ndim == 4:
        model_dict["bcube"] = _normalize_vector_cube_xyzc(np.asarray(model_dict["bcube"]), nx_hint, ny_hint)
    if "chromo_bcube" in model_dict and np.asarray(model_dict["chromo_bcube"]).ndim == 4:
        model_dict["chromo_bcube"] = _normalize_vector_cube_xyzc(
            np.asarray(model_dict["chromo_bcube"]), nx_hint, ny_hint
        )

    for key in ("av_field", "phys_length", "voxel_status", "start_idx", "end_idx"):
        if key in model_dict:
            model_dict[key] = _flatten_status_cube(model_dict[key])

    if "chromo_layers" in model_dict:
        model_dict["chromo_layers"] = int(_scalar_from_any(model_dict["chromo_layers"]))
    if "corona_base" not in model_dict and "bcube" in model_dict and "dz" in model_dict and "chromo_layers" in model_dict:
        bcube_nz = int(np.asarray(model_dict["bcube"]).shape[2])
        full_nz = int(np.asarray(model_dict["dz"]).shape[2])
        chromo_layers = int(model_dict["chromo_layers"])
        inferred = bcube_nz - (full_nz - chromo_layers)
        if 0 <= inferred < bcube_nz:
            model_dict["corona_base"] = inferred
    if "corona_base" in model_dict:
        model_dict["corona_base"] = int(_scalar_from_any(model_dict["corona_base"]))

    required = (
        "dr",
        "dz",
        "bcube",
        "chromo_bcube",
        "chromo_layers",
        "corona_base",
        "chromo_idx",
        "chromo_n",
        "n_p",
        "n_hi",
        "chromo_t",
        "chromo_mask",
        "av_field",
        "phys_length",
        "voxel_status",
        "start_idx",
        "end_idx",
    )
    missing = [key for key in required if key not in model_dict]
    if missing:
        raise KeyError("Missing required CHR datasets after compatibility resolution: " + ", ".join(missing))

    required_header = ("lon", "lat", "obs_time")
    missing_header = [k for k in required_header if k not in header]
    if missing_header:
        raise KeyError("Missing required observation metadata fields: " + ", ".join(missing_header))

    header["obs_time"] = _coerce_time(header["obs_time"])
    return ChromoModelData(header=header, model=model_dict, observer=_normalize_observer_payload(observer_raw))


def _build_chromo_data(file_name):
    """Load chromo model data from H5 or SAV format (auto-detected from filename)."""
    model_payload = load_model(file_name, strict=False)
    return _build_chromo_data_from_payload(model_payload)


def _build_hdf_chromo_data(file_name):
    """Deprecated: Use _build_chromo_data instead. Kept for backwards compatibility."""
    return _build_chromo_data(file_name)


def _build_sav_chromo_data(file_name):
    """Deprecated: Use _build_chromo_data instead. Kept for backwards compatibility."""
    return _build_chromo_data(file_name)


def load_model_dict(model_dict, header):
    lon, lat, obs_time = [header.get(k) for k in ("lon", "lat", "obs_time")]
    obs_time = _coerce_time(obs_time)
    lonc_override = header.get("lonC", None)
    b0sun_override = header.get("b0Sun", None)
    dsun_override = header.get("DSun", None)
    observer_lon_deg = header.get("observer_lon_deg", None)
    observer_lat_deg = header.get("observer_lat_deg", None)
    observer_dsun_cm = header.get("observer_dsun_cm", header.get("DSun", None))
    observer_crln_obs_deg = header.get("observer_crln_obs_deg", header.get("observer_l0_deg", header.get("crln_obs", None)))

    if observer_dsun_cm is not None:
        observer_dsun_cm = float(observer_dsun_cm)
        if observer_dsun_cm < 1e12:
            observer_dsun_cm *= 100.0

    if observer_lon_deg is None or observer_lat_deg is None or observer_dsun_cm is None:
        observer = get_earth(obs_time)
    else:
        observer = SkyCoord(
            lon=float(observer_lon_deg) * u.deg,
            lat=float(observer_lat_deg) * u.deg,
            radius=float(observer_dsun_cm) * u.cm,
            frame=frames.HeliographicStonyhurst,
            obstime=obs_time,
        )

    init_coords = SkyCoord(
        lon * u.deg,
        lat * u.deg,
        frame=frames.HeliographicCarrington,
        rsun=696000 * u.km,
        obstime=obs_time,
        observer=observer,
    )
    hgs = init_coords.transform_to(frames.HeliographicStonyhurst)
    obstime = obs_time.unix - 283996800

    if dsun_override is not None:
        dsun = float(dsun_override)
    else:
        dsun = float(observer_dsun_cm if observer_dsun_cm is not None else sun.earth_distance(obs_time).to(u.cm).value)

    rsun = init_coords.rsun.to(u.cm).value
    if b0sun_override is None:
        if "observer_b0_deg" in header:
            b0sun = float(header["observer_b0_deg"])
        elif "solar_b0" in header:
            b0sun = float(header["solar_b0"])
        elif observer_lat_deg is not None:
            b0sun = float(observer_lat_deg)
        else:
            b0sun = sun.B0(obs_time).value
    else:
        b0sun = float(b0sun_override)

    if lonc_override is None:
        if observer_crln_obs_deg is not None:
            lonc = float(lon) - float(observer_crln_obs_deg)
        else:
            lonc = hgs.lon.to(u.deg).value
    else:
        lonc = float(lonc_override)
    latc = lat

    dr = model_dict["dr"]
    dx = dr[0] * rsun
    dy = dr[1] * rsun
    dz_uniform = dr[2] * rsun
    dz = model_dict["dz"].T * rsun

    s = dz.shape
    sc = model_dict["bcube"].T.shape

    nx = s[2]
    ny = s[1]
    nz = s[0]

    chromo_layers = model_dict["chromo_layers"]
    corona_layers = nz - chromo_layers
    corona_base = model_dict["corona_base"]

    expected_corona_base = int(sc[1] - (nz - chromo_layers))
    if (sc[1] - int(corona_base)) != (nz - chromo_layers):
        if 0 <= expected_corona_base < sc[1]:
            corona_base = expected_corona_base
        else:
            corona_base = int(np.clip(corona_base, 0, max(sc[1] - 1, 0)))

    bx = np.zeros(s, dtype=np.float32, order="C")
    by = np.zeros(s, dtype=np.float32, order="C")
    bz = np.zeros(s, dtype=np.float32, order="C")

    chromo_bcube = model_dict["chromo_bcube"].T
    bcube = model_dict["bcube"].T

    bx[0:chromo_layers, :, :] = chromo_bcube[0, :, :, :]
    by[0:chromo_layers, :, :] = chromo_bcube[1, :, :, :]
    bz[0:chromo_layers, :, :] = chromo_bcube[2, :, :, :]

    cor_target = nz - chromo_layers
    cor_source = max(sc[1] - corona_base, 0)
    cor_copy = min(cor_target, cor_source)
    if cor_copy > 0:
        bx[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[0, corona_base : corona_base + cor_copy, :, :]
        by[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[1, corona_base : corona_base + cor_copy, :, :]
        bz[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[2, corona_base : corona_base + cor_copy, :, :]

    chr_s = (chromo_layers, ny, nx)
    cor_s = (sc[1] - corona_base, ny, nx)
    cor_s2 = (corona_base, ny, nx)

    chromo_n0 = np.zeros(chr_s)
    chromo_np = np.zeros(chr_s)
    chromo_n_hi = np.zeros(chr_s)
    chromo_t0 = np.zeros(chr_s)

    chromo_idx = model_dict["chromo_idx"]
    chromo_n0.flat[chromo_idx] = model_dict["chromo_n"]
    chromo_np.flat[chromo_idx] = model_dict["n_p"]
    chromo_n_hi.flat[chromo_idx] = model_dict["n_hi"]
    chromo_t0.flat[chromo_idx] = model_dict["chromo_t"]

    chromo_mask = model_dict["chromo_mask"].copy()

    voxel_box = {
        "bcube": model_dict["bcube"],
        "dr": model_dict["dr"],
        "chromo_layers": chromo_layers,
        "corona_base": corona_base,
        "chromo_idx": model_dict.get("chromo_idx"),
        "chromo_t": model_dict.get("chromo_t"),
        "chromo_n": model_dict.get("chromo_n"),
    }
    if "start_idx" in model_dict:
        voxel_box["start_idx"] = model_dict["start_idx"]
    voxel_id_xyz = gx_box2id(voxel_box)
    if voxel_id_xyz is None:
        voxel_id = np.zeros(s, dtype=np.uint8)
    else:
        voxel_id = np.asarray(voxel_id_xyz, dtype=np.uint8).transpose((2, 1, 0))
    startidx, endidx = _sanitize_status_mask(model_dict["start_idx"], model_dict["end_idx"], chromo_mask)

    chromo_mask_flat = np.asarray(chromo_mask).ravel(order="C")
    id1 = chromo_mask_flat[startidx]
    id2 = chromo_mask_flat[endidx]

    qb = np.asarray(model_dict["av_field"], dtype=np.float64).reshape(-1)
    ql = np.asarray(model_dict["phys_length"], dtype=np.float64).reshape(-1) * rsun
    uu = np.asarray(model_dict["voxel_status"], dtype=np.uint8).reshape(-1)

    qb[(uu & 4) != 4] = 0
    ql[(uu & 4) != 4] = 0
    id1[(uu & 4) != 4] = 0
    id2[(uu & 4) != 4] = 0

    qb = qb.reshape((nx, ny, sc[1]), order="F")
    ql = ql.reshape((nx, ny, sc[1]), order="F")
    id1 = id1.reshape((nx, ny, sc[1]), order="F")
    id2 = id2.reshape((nx, ny, sc[1]), order="F")

    corona_bavg = qb[:, :, corona_base:sc[1]].T
    chromo_uniform_bavg = qb[:, :, 0:corona_base].T
    corona_l = ql[:, :, corona_base:sc[1]].T
    chromo_uniform_l = ql[:, :, 0:corona_base].T
    corona_id1 = id1[:, :, corona_base:sc[1]].T
    corona_id2 = id2[:, :, corona_base:sc[1]].T
    chromo_uniform_id1 = id1[:, :, 0:corona_base].T
    chromo_uniform_id2 = id2[:, :, 0:corona_base].T

    model_dt_varlist = [
        ("Nx", np.int32),
        ("Ny", np.int32),
        ("Nz", np.int32),
        ("chromo_layers", np.int32),
        ("corona_layers", np.int32),
        ("corona_base", np.int32),
        ("DSun", np.float64),
        ("RSun", np.float64),
        ("b0Sun", np.float64),
        ("lonC", np.float64),
        ("latC", np.float64),
        ("dx", np.float64),
        ("dy", np.float64),
        ("dz_uniform", np.float64),
        ("obstime", np.float64),
        ("dz", np.float32, s),
        ("Bx", np.float32, s),
        ("By", np.float32, s),
        ("Bz", np.float32, s),
        ("chromo_n0", np.float32, chr_s),
        ("chromo_np", np.float32, chr_s),
        ("chromo_nHI", np.float32, chr_s),
        ("chromo_T0", np.float32, chr_s),
        ("corona_Bavg", np.float32, cor_s),
        ("corona_L", np.float32, cor_s),
        ("chromo_uniform_Bavg", np.float32, cor_s2),
        ("chromo_uniform_L", np.float32, cor_s2),
        ("VoxelID", np.byte, s),
        ("corona_ID1", np.byte, cor_s),
        ("corona_ID2", np.byte, cor_s),
        ("chromo_uniform_ID1", np.byte, cor_s2),
        ("chromo_uniform_ID2", np.byte, cor_s2),
    ]

    model_dt = np.dtype(model_dt_varlist)
    model = np.zeros(1, dtype=model_dt)

    locals_map = {
        "Nx": nx,
        "Ny": ny,
        "Nz": nz,
        "chromo_layers": chromo_layers,
        "corona_layers": corona_layers,
        "corona_base": corona_base,
        "DSun": dsun,
        "RSun": rsun,
        "b0Sun": b0sun,
        "lonC": lonc,
        "latC": latc,
        "dx": dx,
        "dy": dy,
        "dz_uniform": dz_uniform,
        "obstime": obstime,
        "dz": dz,
        "Bx": bx,
        "By": by,
        "Bz": bz,
        "chromo_n0": chromo_n0,
        "chromo_np": chromo_np,
        "chromo_nHI": chromo_n_hi,
        "chromo_T0": chromo_t0,
        "corona_Bavg": corona_bavg,
        "corona_L": corona_l,
        "chromo_uniform_Bavg": chromo_uniform_bavg,
        "chromo_uniform_L": chromo_uniform_l,
        "VoxelID": voxel_id,
        "corona_ID1": corona_id1,
        "corona_ID2": corona_id2,
        "chromo_uniform_ID1": chromo_uniform_id1,
        "chromo_uniform_ID2": chromo_uniform_id2,
    }

    for var in model_dt_varlist:
        key, dt = var[0], var[1]
        val = dt(locals_map[key])
        if len(val.shape) > 1:
            model[key] = val.astype(dt, order="C")
        else:
            model[key] = val

    return model, model_dt


def _apply_loader_overrides(
    header,
    DSun=None,
    lonC=None,
    b0Sun=None,
    recompute_observer_ephemeris: bool = False,
    observer_name: str | None = None,
):
    updated = dict(header)
    if recompute_observer_ephemeris:
        recomputed = _recompute_observer_state(updated, _coerce_time(updated["obs_time"]), observer_name or updated.get("observer_name") or updated.get("observer"))
        updated.update(recomputed)
    if DSun is not None:
        updated["DSun"] = float(DSun)
    if lonC is not None:
        updated["lonC"] = float(lonC)
    if b0Sun is not None:
        updated["b0Sun"] = float(b0Sun)
    return updated


def load_model_hdf_with_metadata(
    file_name,
    DSun=None,
    lonC=None,
    b0Sun=None,
    recompute_observer_ephemeris: bool = False,
    observer_name: str | None = None,
):
    data = _build_hdf_chromo_data(file_name)
    file_header = dict(data.header)
    file_header.update(_file_observer_state(file_header))
    header = _apply_loader_overrides(
        file_header,
        DSun=DSun,
        lonC=lonC,
        b0Sun=b0Sun,
        recompute_observer_ephemeris=recompute_observer_ephemeris,
        observer_name=observer_name,
    )
    model, model_dt = load_model_dict(data.model, header)
    header["DSun"] = float(model["DSun"][0])
    header["lonC"] = float(model["lonC"][0])
    header["b0Sun"] = float(model["b0Sun"][0])
    return model, model_dt, dict(header)


def load_model_sav_with_metadata(
    file_name,
    DSun=None,
    lonC=None,
    b0Sun=None,
    recompute_observer_ephemeris: bool = False,
    observer_name: str | None = None,
):
    data = _build_sav_chromo_data(file_name)
    file_header = dict(data.header)
    file_header.update(_file_observer_state(file_header))
    header = _apply_loader_overrides(
        file_header,
        DSun=DSun,
        lonC=lonC,
        b0Sun=b0Sun,
        recompute_observer_ephemeris=recompute_observer_ephemeris,
        observer_name=observer_name,
    )
    model, model_dt = load_model_dict(data.model, header)
    header["DSun"] = float(model["DSun"][0])
    header["lonC"] = float(model["lonC"][0])
    header["b0Sun"] = float(model["b0Sun"][0])
    return model, model_dt, dict(header)


def load_model_with_metadata(
    file_name,
    DSun=None,
    lonC=None,
    b0Sun=None,
    recompute_observer_ephemeris: bool = False,
    observer_name: str | None = None,
):
    data = _build_chromo_data(file_name)
    file_header = dict(data.header)
    file_header.update(_file_observer_state(file_header))
    header = _apply_loader_overrides(
        file_header,
        DSun=DSun,
        lonC=lonC,
        b0Sun=b0Sun,
        recompute_observer_ephemeris=recompute_observer_ephemeris,
        observer_name=observer_name,
    )
    model, model_dt = load_model_dict(data.model, header)
    header["DSun"] = float(model["DSun"][0])
    header["lonC"] = float(model["lonC"][0])
    header["b0Sun"] = float(model["b0Sun"][0])
    return model, model_dt, dict(header), data.observer


def load_model_hdf(file_name, DSun=None, lonC=None, b0Sun=None, recompute_observer_ephemeris: bool = False, observer_name: str | None = None):
    model, model_dt, _header = load_model_hdf_with_metadata(
        file_name,
        DSun=DSun,
        lonC=lonC,
        b0Sun=b0Sun,
        recompute_observer_ephemeris=recompute_observer_ephemeris,
        observer_name=observer_name,
    )
    return model, model_dt


def load_model_sav(file_name, DSun=None, lonC=None, b0Sun=None, recompute_observer_ephemeris: bool = False, observer_name: str | None = None):
    model, model_dt, _header = load_model_sav_with_metadata(
        file_name,
        DSun=DSun,
        lonC=lonC,
        b0Sun=b0Sun,
        recompute_observer_ephemeris=recompute_observer_ephemeris,
        observer_name=observer_name,
    )
    return model, model_dt
