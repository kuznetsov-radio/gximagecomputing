from __future__ import annotations

import re
from pathlib import Path

import astropy.units as u
import h5py
import numpy as np
import scipy.io as io
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.coordinates import frames, get_earth


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
