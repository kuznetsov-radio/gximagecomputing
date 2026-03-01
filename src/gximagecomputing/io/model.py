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
from astropy.time import Time
from sunpy.coordinates import frames, get_earth, sun

from .voxel_id import gx_box2id


def decode_if_bytes(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    return value


def _extract_from_execute(execute_text: str):
    dims_match = re.search(r"--box-dims\s+(\d+)\s+(\d+)\s+(\d+)", execute_text)
    dx_match = re.search(r"--dx-km\s+([0-9.]+)", execute_text)
    box_res_match = re.search(r"--box-res\s+([0-9.]+)", execute_text)
    dims = None
    dx_km = None
    if dims_match:
        dims = (int(dims_match.group(1)), int(dims_match.group(2)), int(dims_match.group(3)))
    if dx_match:
        dx_km = float(dx_match.group(1))
    elif box_res_match:
        dx_km = float(box_res_match.group(1)) * 1000.0
    return dims, dx_km


def extract_center_from_execute(execute_text: str) -> tuple[float, float] | None:
    # Python CLI style: --coords X Y
    m_coords = re.search(
        r"(?:^|\s)--coords(?:\s+|=)([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)",
        execute_text,
    )
    if m_coords:
        return float(m_coords.group(1)), float(m_coords.group(2))

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
    execute_text = None
    if loader_name == "h5":
        with h5py.File(model_path, "r") as f:
            if "metadata" in f and "execute" in f["metadata"]:
                execute_text = decode_if_bytes(f["metadata"]["execute"][()])
    else:
        data = io.readsav(str(model_path))
        box = data.box
        if "EXECUTE" in box.dtype.names:
            execute_text = decode_if_bytes(box.execute[0])

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
    execute_text = None
    if loader_name == "h5":
        with h5py.File(model_path, "r") as f:
            if "metadata" in f and "execute" in f["metadata"]:
                execute_text = decode_if_bytes(f["metadata"]["execute"][()])
    else:
        data = io.readsav(str(model_path))
        box = data.box
        if "EXECUTE" in box.dtype.names:
            execute_text = decode_if_bytes(box.execute[0])

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


def estimate_hpc_center(model) -> tuple[float, float]:
    unix = float(model["obstime"][0]) + 283996800.0
    obs_time = Time(unix, format="unix")
    observer = get_earth(obs_time)
    center_hgs = SkyCoord(
        lon=float(model["lonC"][0]) * u.deg,
        lat=float(model["latC"][0]) * u.deg,
        radius=(float(model["RSun"][0]) / 1e5) * u.km,
        frame=frames.HeliographicStonyhurst,
        obstime=obs_time,
        observer=observer,
    )
    center_hpc = center_hgs.transform_to(frames.Helioprojective(observer=observer, obstime=obs_time))
    return float(center_hpc.Tx.to_value(u.arcsec)), float(center_hpc.Ty.to_value(u.arcsec))


@dataclass
class ChromoModelData:
    header: Dict[str, Any]
    model: Dict[str, np.ndarray]


def _coerce_time(obs_time):
    if isinstance(obs_time, Time):
        return obs_time
    if isinstance(obs_time, (bytes, np.bytes_)):
        obs_time = obs_time.decode("utf-8")
    return Time(obs_time)


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


def _fill_header_from_base_index(model_f: h5py.File, header: Dict[str, Any]) -> None:
    if "base" not in model_f or "index" not in model_f["base"]:
        return
    raw = _decode_dataset_scalar(model_f["base"]["index"])
    text = raw if isinstance(raw, str) else str(raw)

    def _extract_fits_card(key: str) -> str | None:
        m = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([^\n/]+)", text)
        if not m:
            return None
        value = m.group(1).strip()
        if value.startswith("'"):
            m_str = re.match(r"'([^']*)'", value)
            if m_str:
                return m_str.group(1).strip()
        return value.strip().strip("'")

    if "lon" not in header:
        mlon = re.search(r",\s*([+-]?\d+(?:\.\d+)?)\s*,\s*b'CRLN-CEA'", text)
        if mlon:
            header["lon"] = float(mlon.group(1))
        else:
            crval1 = _extract_fits_card("CRVAL1")
            if crval1 is not None:
                try:
                    header["lon"] = float(crval1)
                except ValueError:
                    pass
    if "lat" not in header:
        mlat = re.search(r",\s*([+-]?\d+(?:\.\d+)?)\s*,\s*b'CRLT-CEA'", text)
        if mlat:
            header["lat"] = float(mlat.group(1))
        else:
            crval2 = _extract_fits_card("CRVAL2")
            if crval2 is not None:
                try:
                    header["lat"] = float(crval2)
                except ValueError:
                    pass
    if "obs_time" not in header:
        date_obs = _extract_fits_card("DATE-OBS") or _extract_fits_card("DATE_OBS")
        if date_obs:
            header["obs_time"] = date_obs
        else:
            mtime = re.search(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?", text)
            if mtime:
                header["obs_time"] = mtime.group(0)
    if "dsun_obs" not in header:
        dsun_obs = _extract_fits_card("DSUN_OBS")
        if dsun_obs is not None:
            try:
                header["dsun_obs"] = float(dsun_obs)
            except ValueError:
                pass
        else:
            mds = re.search(
                r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?'\s*,\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)",
                text,
            )
            if mds:
                header["dsun_obs"] = float(mds.group(1))
    # Do not infer lonC from index fallback tags (HGLN_OBS/CRVAL1).
    # Keep lonC derived from coordinate transforms in load_model_dict,
    # unless an explicit lonC override is provided by higher-level callers.


def _build_hdf_chromo_data(file_name):
    with h5py.File(file_name, "r") as model_f:
        if "chromo" not in model_f:
            raise KeyError("Missing required '/chromo' group in HDF5 model.")
        chromo_box = model_f["chromo"]
        model_dict = _load_h5_group(chromo_box)
        header = {k: decode_if_bytes(v) for k, v in dict(chromo_box.attrs).items()}

        ny_hint = None
        nx_hint = None
        if "base" in model_f and "bx" in model_f["base"]:
            base_shape = np.asarray(model_f["base"]["bx"][:]).shape
            if len(base_shape) == 2:
                ny_hint, nx_hint = int(base_shape[0]), int(base_shape[1])

        if "lines" in model_f:
            lines_box = _load_h5_group(model_f["lines"])
            for key in ("av_field", "phys_length", "voxel_status", "start_idx", "end_idx"):
                if key not in model_dict and key in lines_box:
                    model_dict[key] = lines_box[key]

        if "corona" in model_f:
            corona_box = _load_h5_group(model_f["corona"])
            if "dr" not in model_dict and "dr" in corona_box:
                model_dict["dr"] = corona_box["dr"]
            if "corona_base" not in model_dict and "corona_base" in corona_box:
                model_dict["corona_base"] = corona_box["corona_base"]
            if "bcube" not in model_dict:
                if all(k in corona_box for k in ("bx", "by", "bz")):
                    model_dict["bcube"] = _components_to_vector_cube_xyzc(
                        corona_box["bx"], corona_box["by"], corona_box["bz"], nx_hint, ny_hint
                    )
                elif "bcube" in corona_box:
                    model_dict["bcube"] = corona_box["bcube"]

        if "dr" not in model_dict and "dr" in chromo_box:
            model_dict["dr"] = chromo_box["dr"][:]
        if "corona_base" not in model_dict and "corona_base" in chromo_box:
            model_dict["corona_base"] = chromo_box["corona_base"][()]
        if "chromo_bcube" not in model_dict:
            if all(k in model_dict for k in ("bx", "by", "bz")):
                model_dict["chromo_bcube"] = _components_to_vector_cube_xyzc(
                    model_dict["bx"], model_dict["by"], model_dict["bz"], nx_hint, ny_hint
                )
            elif "chromo_bcube" in chromo_box:
                model_dict["chromo_bcube"] = chromo_box["chromo_bcube"][:]
        if "bcube" not in model_dict and "bcube" in chromo_box:
            model_dict["bcube"] = chromo_box["bcube"][:]

        if "chromo_mask" not in model_dict and "base" in model_f and "chromo_mask" in model_f["base"]:
            model_dict["chromo_mask"] = model_f["base"]["chromo_mask"][:]
        # Keep lonC computation centralized in load_model_dict unless explicitly overridden.
        if "corona" in model_f:
            for key, value in dict(model_f["corona"].attrs).items():
                if key not in header:
                    header[key] = decode_if_bytes(value)
        if "metadata" in model_f and "execute" in model_f["metadata"] and "obs_time" not in header:
            execute_text = _decode_dataset_scalar(model_f["metadata"]["execute"])
            if isinstance(execute_text, str):
                t_fallback = _parse_execute_time(execute_text)
                if t_fallback:
                    header["obs_time"] = t_fallback
        _fill_header_from_base_index(model_f, header)

    if "dz" in model_dict and np.asarray(model_dict["dz"]).ndim == 3:
        model_dict["dz"] = _normalize_cube_xyz(np.asarray(model_dict["dz"]), nx_hint, ny_hint)
    if "bcube" in model_dict and np.asarray(model_dict["bcube"]).ndim == 4:
        model_dict["bcube"] = _normalize_vector_cube_xyzc(np.asarray(model_dict["bcube"]), nx_hint, ny_hint)
    if "chromo_bcube" in model_dict and np.asarray(model_dict["chromo_bcube"]).ndim == 4:
        model_dict["chromo_bcube"] = _normalize_vector_cube_xyzc(
            np.asarray(model_dict["chromo_bcube"]), nx_hint, ny_hint
        )
    if "chromo_layers" in model_dict:
        model_dict["chromo_layers"] = int(_scalar_from_any(model_dict["chromo_layers"]))
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
    return ChromoModelData(header=header, model=model_dict)


def _build_sav_chromo_data(file_name):
    model_data = io.readsav(file_name)
    box = model_data.box
    lon = box.index[0].CRVAL1[0]
    lat = box.index[0].CRVAL2[0]
    aptime = Time(box.index[0]["DATE_OBS"][0])
    header = {
        "lon": lon,
        "lat": lat,
        "dsun_obs": box.index[0]["DSUN_OBS"][0],
        "obs_time": aptime,
    }
    # Do not seed lonC directly from INDEX/HGLN_OBS here; keep model lonC
    # derived in load_model_dict unless explicitly provided by caller.

    model_dict = {
        "dr": box.dr[0],
        "dz": box.dz[0].transpose((2, 1, 0)),
        "bcube": box.bcube[0].transpose((3, 2, 1, 0)),
        "chromo_bcube": box.chromo_bcube[0].transpose((3, 2, 1, 0)),
        "chromo_layers": box.chromo_layers[0],
        "corona_base": box.corona_base[0],
        "chromo_idx": box.chromo_idx[0].astype(np.int64),
        "chromo_n": box.chromo_n[0],
        "n_p": box.n_p[0],
        "n_hi": box.n_hi[0],
        "chromo_t": box.chromo_t[0],
        "chromo_mask": box["base"][0]["chromo_mask"][0],
    }

    if "AVFIELD" in box.dtype.names:
        model_dict["av_field"] = box.avfield[0].T.reshape(-1, order="F")
        model_dict["phys_length"] = box.physlength[0].T.reshape(-1, order="F")
        model_dict["voxel_status"] = box.status[0].T.reshape(-1, order="F")
        model_dict["start_idx"] = box.startidx[0].T.reshape(-1, order="F")
        model_dict["end_idx"] = box.endidx[0].T.reshape(-1, order="F")
    else:
        sc = box.bcube[0].shape
        nx, ny = sc[3], sc[2]
        qb = np.zeros((nx, ny, sc[1]), dtype=np.float64)
        ql = np.zeros((nx, ny, sc[1]), dtype=np.float64)
        uu = np.zeros((nx, ny, sc[1]), dtype=np.uint8)
        idx = np.unravel_index(box.idx[0], qb.shape, order="F")
        qb[idx] = box.bmed[0]
        ql[idx] = box.length[0]
        uu[idx] = 4
        model_dict["av_field"] = qb.reshape(-1, order="F")
        model_dict["phys_length"] = ql.reshape(-1, order="F")
        model_dict["voxel_status"] = uu.reshape(-1, order="F")
        model_dict["start_idx"] = np.zeros(qb.size, dtype=np.int64)
        model_dict["end_idx"] = np.zeros(qb.size, dtype=np.int64)

    return ChromoModelData(header=header, model=model_dict)


def load_model_dict(model_dict, header):
    lon, lat, obs_time = [header.get(k) for k in ("lon", "lat", "obs_time")]
    dsun_obs = header.get("dsun_obs", None)
    obs_time = _coerce_time(obs_time)
    lonc_override = header.get("lonC", None)
    b0sun_override = header.get("b0Sun", None)
    dsun_override = header.get("DSun", None)
    init_coords = SkyCoord(
        lon * u.deg, lat * u.deg, frame=frames.HeliographicCarrington, rsun=696000 * u.km, obstime=obs_time, observer="earth"
    )
    hgs = init_coords.transform_to(frames.HeliographicStonyhurst)
    obstime = obs_time.unix - 283996800

    if dsun_override is not None:
        dsun = float(dsun_override)
    else:
        # Match IDL get_sun(anytim(obsTime,/ex)) behavior for the default Earth observer.
        dsun = sun.earth_distance(obs_time).to(u.cm).value
        if (not np.isfinite(dsun)) and (dsun_obs is not None):
            # Last-resort fallback for malformed/unsupported times.
            dsun = float(dsun_obs) * 1e2
    # Match updated IDL LoadGXmodel: normalize DSUN to centimeters when header/HDF stores meters.
    if dsun < 1e12:
        dsun *= 100.0

    rsun = init_coords.rsun.to(u.cm).value
    b0sun = sun.B0(obs_time).value if b0sun_override is None else float(b0sun_override)

    lonc = hgs.lon.to(u.deg).value if lonc_override is None else float(lonc_override)
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


def load_model_hdf(file_name):
    data = _build_hdf_chromo_data(file_name)
    return load_model_dict(data.model, data.header)


def load_model_sav(file_name):
    data = _build_sav_chromo_data(file_name)
    return load_model_dict(data.model, data.header)
