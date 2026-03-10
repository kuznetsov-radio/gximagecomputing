from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from astropy.io import fits


def save_h5_maps(
    result: dict,
    freqlist: list[float],
    out_h5: Path,
    model_path: Path,
    model_format: str,
    xc: float,
    yc: float,
    dx: float,
    dy: float,
    obs_time_iso: str,
) -> Path:
    ti = np.asarray(result["TI"], dtype=np.float32)  # (ny, nx, nf)
    tv = np.asarray(result["TV"], dtype=np.float32)  # (ny, nx, nf)
    ny, nx, nf = ti.shape

    # Requested container layout: [nx, ny, nf, 2] where last axis is [TI, TV].
    ti_nxnyf = np.transpose(ti, (1, 0, 2))
    tv_nxnyf = np.transpose(tv, (1, 0, 2))
    cube = np.stack([ti_nxnyf, tv_nxnyf], axis=-1)  # (nx, ny, nf, 2)

    hdr = fits.Header()
    hdr["SIMPLE"] = 1
    hdr["BITPIX"] = -32
    hdr["NAXIS"] = 4
    hdr["NAXIS1"] = int(nx)
    hdr["NAXIS2"] = int(ny)
    hdr["NAXIS3"] = int(nf)
    hdr["NAXIS4"] = 2
    hdr["CTYPE1"] = "HPLN-TAN"
    hdr["CTYPE2"] = "HPLT-TAN"
    hdr["CUNIT1"] = "arcsec"
    hdr["CUNIT2"] = "arcsec"
    hdr["CDELT1"] = float(dx)
    hdr["CDELT2"] = float(dy)
    hdr["CRPIX1"] = (nx + 1.0) / 2.0
    hdr["CRPIX2"] = (ny + 1.0) / 2.0
    hdr["CRVAL1"] = float(xc)
    hdr["CRVAL2"] = float(yc)
    hdr["CTYPE3"] = "FREQ"
    hdr["CUNIT3"] = "GHz"
    hdr["CRPIX3"] = 1.0
    hdr["CRVAL3"] = float(freqlist[0]) if len(freqlist) else 0.0
    hdr["CDELT3"] = float(freqlist[1] - freqlist[0]) if len(freqlist) > 1 else 0.0
    hdr["CTYPE4"] = "STOKES"
    hdr["CUNIT4"] = ""
    hdr["CRPIX4"] = 1.0
    hdr["CRVAL4"] = 1.0
    hdr["CDELT4"] = 1.0
    hdr["DATE-OBS"] = obs_time_iso
    hdr["BUNIT"] = "K"

    map_ids = [f"Tb_I {f:.2f} GHz" for f in freqlist] + [f"Tb_V {f:.2f} GHz" for f in freqlist]
    map_ids_arr = np.asarray(map_ids, dtype="S64")
    stokes_ids_arr = np.asarray(["TI", "TV"], dtype="S8")

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        maps = f.create_group("maps")
        maps.create_dataset("data", data=cube, compression="gzip", compression_opts=4)
        maps.create_dataset("freqlist_ghz", data=np.asarray(freqlist, dtype=np.float64))
        maps.create_dataset("stokes_ids", data=stokes_ids_arr)
        maps.create_dataset("map_ids", data=map_ids_arr)

        meta = f.create_group("metadata")
        meta.create_dataset("index_header", data=np.bytes_(hdr.tostring(sep="\n", endcard=True)))
        meta.create_dataset("model_path", data=np.bytes_(str(model_path)))
        meta.create_dataset("model_format", data=np.bytes_(model_format))
        meta.create_dataset("date_obs", data=np.bytes_(obs_time_iso))
        meta.create_dataset("xc_arcsec", data=float(xc))
        meta.create_dataset("yc_arcsec", data=float(yc))
        meta.create_dataset("dx_arcsec", data=float(dx))
        meta.create_dataset("dy_arcsec", data=float(dy))
        meta.create_dataset("nx", data=int(nx))
        meta.create_dataset("ny", data=int(ny))
        meta.create_dataset("nf", data=int(nf))
    return out_h5


def save_h5_euv_maps(
    flux_corona: np.ndarray,
    flux_tr: np.ndarray,
    channels: list[str],
    out_h5: Path,
    model_path: Path,
    model_format: str,
    xc: float,
    yc: float,
    dx: float,
    dy: float,
    obs_time_iso: str,
    instrument: str = "AIA",
) -> Path:
    cor = np.asarray(flux_corona, dtype=np.float32)  # (ny, nx, nch)
    tr = np.asarray(flux_tr, dtype=np.float32)  # (ny, nx, nch)
    if cor.shape != tr.shape:
        raise ValueError(f"EUV corona/TR shape mismatch: {cor.shape} vs {tr.shape}")
    ny, nx, nch = cor.shape
    if len(channels) != nch:
        raise ValueError(f"Channel count mismatch: channels={len(channels)} vs cube nch={nch}")

    cor_nxnych = np.transpose(cor, (1, 0, 2))
    tr_nxnych = np.transpose(tr, (1, 0, 2))
    cube = np.stack([cor_nxnych, tr_nxnych], axis=-1)  # (nx, ny, nch, 2)

    hdr = fits.Header()
    hdr["SIMPLE"] = 1
    hdr["BITPIX"] = -32
    hdr["NAXIS"] = 4
    hdr["NAXIS1"] = int(nx)
    hdr["NAXIS2"] = int(ny)
    hdr["NAXIS3"] = int(nch)
    hdr["NAXIS4"] = 2
    hdr["CTYPE1"] = "HPLN-TAN"
    hdr["CTYPE2"] = "HPLT-TAN"
    hdr["CUNIT1"] = "arcsec"
    hdr["CUNIT2"] = "arcsec"
    hdr["CDELT1"] = float(dx)
    hdr["CDELT2"] = float(dy)
    hdr["CRPIX1"] = (nx + 1.0) / 2.0
    hdr["CRPIX2"] = (ny + 1.0) / 2.0
    hdr["CRVAL1"] = float(xc)
    hdr["CRVAL2"] = float(yc)
    hdr["CTYPE3"] = "CHANNEL"
    hdr["CUNIT3"] = ""
    hdr["CRPIX3"] = 1.0
    hdr["CRVAL3"] = 1.0
    hdr["CDELT3"] = 1.0
    hdr["CTYPE4"] = "COMPONENT"
    hdr["CUNIT4"] = ""
    hdr["CRPIX4"] = 1.0
    hdr["CRVAL4"] = 1.0
    hdr["CDELT4"] = 1.0
    hdr["DATE-OBS"] = obs_time_iso
    hdr["BUNIT"] = "DN s^-1 pix^-1"

    map_ids = [f"{instrument} {c} COR" for c in channels] + [f"{instrument} {c} TR" for c in channels]
    comp_ids = np.asarray(["CORONA", "TR"], dtype="S16")
    ch_ids = np.asarray([str(c) for c in channels], dtype="S16")
    map_ids_arr = np.asarray(map_ids, dtype="S64")

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        maps = f.create_group("maps")
        maps.create_dataset("data", data=cube, compression="gzip", compression_opts=4)
        maps.create_dataset("channel_ids", data=ch_ids)
        maps.create_dataset("component_ids", data=comp_ids)
        maps.create_dataset("map_ids", data=map_ids_arr)

        meta = f.create_group("metadata")
        meta.create_dataset("index_header", data=np.bytes_(hdr.tostring(sep="\n", endcard=True)))
        meta.create_dataset("model_path", data=np.bytes_(str(model_path)))
        meta.create_dataset("model_format", data=np.bytes_(model_format))
        meta.create_dataset("instrument", data=np.bytes_(str(instrument)))
        meta.create_dataset("date_obs", data=np.bytes_(obs_time_iso))
        meta.create_dataset("xc_arcsec", data=float(xc))
        meta.create_dataset("yc_arcsec", data=float(yc))
        meta.create_dataset("dx_arcsec", data=float(dx))
        meta.create_dataset("dy_arcsec", data=float(dy))
        meta.create_dataset("nx", data=int(nx))
        meta.create_dataset("ny", data=int(ny))
        meta.create_dataset("nch", data=int(nch))
    return out_h5
