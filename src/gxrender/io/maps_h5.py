from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from astropy.io import fits


def _header_text(header: fits.Header | dict) -> str:
    hdr = header if isinstance(header, fits.Header) else fits.Header(header)
    return hdr.tostring(sep="\n", endcard=True)


def _write_common_metadata(
    meta: h5py.Group,
    *,
    wcs_header: fits.Header | dict | None,
    observer_name: str | None,
    observer_source: str | None,
    observer_warnings: list[str] | tuple[str, ...] | None,
    l0_deg: float | None,
    b0_deg: float | None,
    dsun_cm: float | None,
    rsun_cm: float | None,
    rsun_arcsec: float | None,
) -> None:
    if wcs_header is not None:
        header_text = _header_text(wcs_header)
        meta.create_dataset("wcs_header", data=np.bytes_(header_text))
        meta.create_dataset("index_header", data=np.bytes_(header_text))
    if observer_name is not None:
        meta.create_dataset("observer_name", data=np.bytes_(str(observer_name)))
    if observer_source is not None:
        meta.create_dataset("observer_source", data=np.bytes_(str(observer_source)))
    if observer_warnings:
        meta.create_dataset("observer_warnings", data=np.bytes_("\n".join(str(v) for v in observer_warnings)))
    if l0_deg is not None:
        meta.create_dataset("observer_l0_deg", data=float(l0_deg))
    if b0_deg is not None:
        meta.create_dataset("observer_b0_deg", data=float(b0_deg))
    if dsun_cm is not None:
        meta.create_dataset("observer_dsun_cm", data=float(dsun_cm))
    if rsun_cm is not None:
        meta.create_dataset("observer_rsun_cm", data=float(rsun_cm))
    if rsun_arcsec is not None:
        meta.create_dataset("observer_rsun_arcsec", data=float(rsun_arcsec))


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
    wcs_header: fits.Header | dict | None = None,
    observer_name: str | None = None,
    observer_source: str | None = None,
    observer_warnings: list[str] | tuple[str, ...] | None = None,
    l0_deg: float | None = None,
    b0_deg: float | None = None,
    dsun_cm: float | None = None,
    rsun_cm: float | None = None,
    rsun_arcsec: float | None = None,
    tbase: float | None = None,
    nbase: float | None = None,
    q0: float | None = None,
    a: float | None = None,
    b: float | None = None,
    corona_mode: int | None = None,
    shtable: np.ndarray | None = None,
) -> Path:
    ti = np.asarray(result["TI"], dtype=np.float32)  # (ny, nx, nf)
    tv = np.asarray(result["TV"], dtype=np.float32)  # (ny, nx, nf)
    ny, nx, nf = ti.shape

    # Requested container layout: [nx, ny, nf, 2] where last axis is [TI, TV].
    ti_nxnyf = np.transpose(ti, (1, 0, 2))
    tv_nxnyf = np.transpose(tv, (1, 0, 2))
    cube = np.stack([ti_nxnyf, tv_nxnyf], axis=-1)  # (nx, ny, nf, 2)

    map_ids = [f"Tb_I {f:.2f} GHz" for f in freqlist] + [f"Tb_V {f:.2f} GHz" for f in freqlist]
    map_ids_arr = np.asarray(map_ids, dtype="S64")
    stokes_ids_arr = np.asarray(["TI", "TV"], dtype="S8")
    if wcs_header is None:
        wcs_header = fits.Header(
            {
                "NAXIS": 2,
                "NAXIS1": int(nx),
                "NAXIS2": int(ny),
                "CTYPE1": "HPLN-TAN",
                "CTYPE2": "HPLT-TAN",
                "CUNIT1": "arcsec",
                "CUNIT2": "arcsec",
                "CDELT1": float(dx),
                "CDELT2": float(dy),
                "CRPIX1": (nx + 1.0) / 2.0,
                "CRPIX2": (ny + 1.0) / 2.0,
                "CRVAL1": float(xc),
                "CRVAL2": float(yc),
                "DATE-OBS": obs_time_iso,
                "BUNIT": "K",
            }
        )

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        maps = f.create_group("maps")
        maps.create_dataset("data", data=cube, compression="gzip", compression_opts=4)
        maps.create_dataset("freqlist_ghz", data=np.asarray(freqlist, dtype=np.float64))
        maps.create_dataset("stokes_ids", data=stokes_ids_arr)
        maps.create_dataset("map_ids", data=map_ids_arr)

        meta = f.create_group("metadata")
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
        if tbase is not None:
            meta.create_dataset("tbase_k", data=float(tbase))
        if nbase is not None:
            meta.create_dataset("nbase_cm3", data=float(nbase))
        if q0 is not None:
            meta.create_dataset("q0", data=float(q0))
        if a is not None:
            meta.create_dataset("a", data=float(a))
        if b is not None:
            meta.create_dataset("b", data=float(b))
        if corona_mode is not None:
            meta.create_dataset("corona_mode", data=int(corona_mode))
        if shtable is not None:
            meta.create_dataset("shtable", data=np.asarray(shtable, dtype=np.float64))
        _write_common_metadata(
            meta,
            wcs_header=wcs_header,
            observer_name=observer_name,
            observer_source=observer_source,
            observer_warnings=observer_warnings,
            l0_deg=l0_deg,
            b0_deg=b0_deg,
            dsun_cm=dsun_cm,
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )
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
    wcs_header: fits.Header | dict | None = None,
    observer_name: str | None = None,
    observer_source: str | None = None,
    observer_warnings: list[str] | tuple[str, ...] | None = None,
    l0_deg: float | None = None,
    b0_deg: float | None = None,
    dsun_cm: float | None = None,
    rsun_cm: float | None = None,
    rsun_arcsec: float | None = None,
    tbase: float | None = None,
    nbase: float | None = None,
    q0: float | None = None,
    a: float | None = None,
    b: float | None = None,
    corona_mode: int | None = None,
    shtable: np.ndarray | None = None,
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

    map_ids = [f"{instrument} {c} COR" for c in channels] + [f"{instrument} {c} TR" for c in channels]
    comp_ids = np.asarray(["CORONA", "TR"], dtype="S16")
    ch_ids = np.asarray([str(c) for c in channels], dtype="S16")
    map_ids_arr = np.asarray(map_ids, dtype="S64")
    if wcs_header is None:
        wcs_header = fits.Header(
            {
                "NAXIS": 2,
                "NAXIS1": int(nx),
                "NAXIS2": int(ny),
                "CTYPE1": "HPLN-TAN",
                "CTYPE2": "HPLT-TAN",
                "CUNIT1": "arcsec",
                "CUNIT2": "arcsec",
                "CDELT1": float(dx),
                "CDELT2": float(dy),
                "CRPIX1": (nx + 1.0) / 2.0,
                "CRPIX2": (ny + 1.0) / 2.0,
                "CRVAL1": float(xc),
                "CRVAL2": float(yc),
                "DATE-OBS": obs_time_iso,
                "BUNIT": "DN s^-1 pix^-1",
            }
        )

    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        maps = f.create_group("maps")
        maps.create_dataset("data", data=cube, compression="gzip", compression_opts=4)
        maps.create_dataset("channel_ids", data=ch_ids)
        maps.create_dataset("component_ids", data=comp_ids)
        maps.create_dataset("map_ids", data=map_ids_arr)

        meta = f.create_group("metadata")
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
        if tbase is not None:
            meta.create_dataset("tbase_k", data=float(tbase))
        if nbase is not None:
            meta.create_dataset("nbase_cm3", data=float(nbase))
        if q0 is not None:
            meta.create_dataset("q0", data=float(q0))
        if a is not None:
            meta.create_dataset("a", data=float(a))
        if b is not None:
            meta.create_dataset("b", data=float(b))
        if corona_mode is not None:
            meta.create_dataset("corona_mode", data=int(corona_mode))
        if shtable is not None:
            meta.create_dataset("shtable", data=np.asarray(shtable, dtype=np.float64))
        _write_common_metadata(
            meta,
            wcs_header=wcs_header,
            observer_name=observer_name,
            observer_source=observer_source,
            observer_warnings=observer_warnings,
            l0_deg=l0_deg,
            b0_deg=b0_deg,
            dsun_cm=dsun_cm,
            rsun_cm=rsun_cm,
            rsun_arcsec=rsun_arcsec,
        )
    return out_h5
