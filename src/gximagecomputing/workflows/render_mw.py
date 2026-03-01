from __future__ import annotations

import argparse
from pathlib import Path

import gximagecomputing
import h5py
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.io import fits

from gximagecomputing.io.maps_h5 import save_h5_maps
from gximagecomputing.workflows._render_common import DEFAULT_OUTDIR, model_obstime_iso, plasma_defaults, prepare_common_inputs


def _build_2d_wcs_meta(index_header: str, nx: int, ny: int, date_obs: str) -> dict:
    meta = {
        "naxis": 2,
        "naxis1": int(nx),
        "naxis2": int(ny),
        "ctype1": "HPLN-TAN",
        "ctype2": "HPLT-TAN",
        "cunit1": "arcsec",
        "cunit2": "arcsec",
        "cdelt1": 1.0,
        "cdelt2": 1.0,
        "crpix1": (nx + 1.0) / 2.0,
        "crpix2": (ny + 1.0) / 2.0,
        "crval1": 0.0,
        "crval2": 0.0,
        "date-obs": date_obs or "",
        "bunit": "K",
    }
    if not index_header:
        return meta
    try:
        hdr = fits.Header.fromstring(index_header, sep="\n")
    except Exception:
        return meta
    for k in ("CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CDELT1", "CDELT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "BUNIT"):
        if k in hdr:
            meta[k.lower()] = hdr[k]
    if "DATE-OBS" in hdr and not meta["date-obs"]:
        meta["date-obs"] = hdr["DATE-OBS"]
    return meta


def save_preview(h5_path: Path, out_png: Path, title: str) -> None:
    with h5py.File(h5_path, "r") as f:
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # [nx, ny, nf, 2]
        freqs = np.asarray(f["maps"]["freqlist_ghz"], dtype=np.float64)
        index_header = ""
        date_obs = ""
        if "metadata" in f and "index_header" in f["metadata"]:
            raw = f["metadata"]["index_header"][()]
            index_header = raw.decode("utf-8", errors="replace") if isinstance(raw, (bytes, np.bytes_)) else str(raw)
        if "metadata" in f and "date_obs" in f["metadata"]:
            raw_date = f["metadata"]["date_obs"][()]
            date_obs = raw_date.decode("utf-8", errors="replace") if isinstance(raw_date, (bytes, np.bytes_)) else str(raw_date)

    if cube.ndim != 4 or cube.shape[-1] != 2:
        raise ValueError(f"Unexpected rendered-map cube shape: {cube.shape}")

    # maps/data is [nx, ny, nf, stokes]; SunPy map data expects [ny, nx].
    ti = cube[:, :, 0, 0].T
    tv = cube[:, :, 0, 1].T
    ny, nx = ti.shape
    meta_base = _build_2d_wcs_meta(index_header=index_header, nx=nx, ny=ny, date_obs=date_obs)
    freq0 = float(freqs[0]) if len(freqs) else float("nan")

    ti_meta = dict(meta_base)
    ti_meta["content"] = "TI"
    ti_meta["freqghz"] = freq0
    tv_meta = dict(meta_base)
    tv_meta["content"] = "TV"
    tv_meta["freqghz"] = freq0

    m_ti = sunpy.map.Map(ti, ti_meta)
    m_tv = sunpy.map.Map(tv, tv_meta)

    fig = plt.figure(figsize=(10, 4))
    ax_ti = fig.add_subplot(1, 2, 1, projection=m_ti)
    ax_tv = fig.add_subplot(1, 2, 2, projection=m_tv)
    im_ti = m_ti.plot(axes=ax_ti, cmap="inferno", interpolation="nearest")
    im_tv = m_tv.plot(axes=ax_tv, cmap="coolwarm", interpolation="nearest")
    ax_ti.set_title(f"{title} TI @ {freq0:.2f} GHz\\n{date_obs}")
    ax_tv.set_title(f"{title} TV @ {freq0:.2f} GHz\\n{date_obs}")
    ax_ti.set_xlabel("Solar X [arcsec]")
    ax_ti.set_ylabel("Solar Y [arcsec]")
    ax_tv.set_xlabel("Solar X [arcsec]")
    ax_tv.set_ylabel("Solar Y [arcsec]")
    cbar0 = fig.colorbar(im_ti, ax=ax_ti, orientation="vertical", fraction=0.046, pad=0.04)
    cbar1 = fig.colorbar(im_tv, ax=ax_tv, orientation="vertical", fraction=0.046, pad=0.04)
    cbar0.set_label("Brightness Temperature [K]")
    cbar1.set_label("Brightness Temperature [K]")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def save_sunpy_maps(
    result: dict,
    freqlist: list[float],
    out_dir: Path,
    stem: str,
    xc: float,
    yc: float,
    dx: float,
    dy: float,
    obs_time_iso: str,
) -> list[str]:
    ti = np.asarray(result["TI"], dtype=np.float32)
    tv = np.asarray(result["TV"], dtype=np.float32)
    ny, nx, nf = ti.shape
    files: list[str] = []
    for i in range(nf):
        freq = float(freqlist[i])
        base_meta = {
            "naxis": 2,
            "naxis1": int(nx),
            "naxis2": int(ny),
            "ctype1": "HPLN-TAN",
            "ctype2": "HPLT-TAN",
            "cunit1": "arcsec",
            "cunit2": "arcsec",
            "cdelt1": float(dx),
            "cdelt2": float(dy),
            "crpix1": (nx + 1.0) / 2.0,
            "crpix2": (ny + 1.0) / 2.0,
            "crval1": float(xc),
            "crval2": float(yc),
            "date-obs": obs_time_iso,
            "freqghz": freq,
            "bunit": "K",
        }
        meta_ti = dict(base_meta)
        meta_ti["content"] = "TI"
        meta_tv = dict(base_meta)
        meta_tv["content"] = "TV"

        m_ti = sunpy.map.Map(ti[:, :, i], meta_ti)
        m_tv = sunpy.map.Map(tv[:, :, i], meta_tv)
        ti_path = out_dir / f"{stem}.ti.{freq:05.1f}GHz.fits"
        tv_path = out_dir / f"{stem}.tv.{freq:05.1f}GHz.fits"
        m_ti.save(str(ti_path), overwrite=True)
        m_tv.save(str(tv_path), overwrite=True)
        files.extend([str(ti_path), str(tv_path)])
    return files


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render MW maps from a single H5/SAV CHR model.")
    p.add_argument("--model-path", type=Path, required=True, help="Path to input CHR model (.h5 or .sav).")
    p.add_argument("--model-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument(
        "--ebtel-path",
        type=str,
        default=None,
        help='Optional EBTEL table (.sav). Use "" to explicitly disable and force isothermal mode.',
    )
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument(
        "--output-name",
        type=str,
        default=None,
        help="Default: <model-basename>_py_mw_maps.h5",
    )
    p.add_argument(
        "--output-format",
        choices=["h5", "fits", "both"],
        default="h5",
        help="h5: compact cube container (default); fits: per-frequency FITS maps; both: h5+fits.",
    )
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--dsun-cm", type=float, default=None, help="Override model.DSun before rendering (cm).")
    p.add_argument("--lonc-deg", type=float, default=None, help="Override model.lonC before rendering (deg).")
    p.add_argument("--b0sun-deg", type=float, default=None, help="Override model.b0Sun before rendering (deg).")
    p.add_argument("--dx", type=float, default=None, help="Pixel scale along X (arcsec/pixel). Default: 2.0")
    p.add_argument("--dy", type=float, default=None, help="Pixel scale along Y (arcsec/pixel). Default: 2.0")
    p.add_argument(
        "--pixel-scale-arcsec",
        type=float,
        default=None,
        help="Legacy shorthand: sets both --dx and --dy when they are omitted.",
    )
    p.add_argument("--nx", type=int, default=None, help="Optional fixed output width in pixels.")
    p.add_argument("--ny", type=int, default=None, help="Optional fixed output height in pixels.")
    p.add_argument("--xrange", type=float, nargs=2, default=None, metavar=("XMIN", "XMAX"))
    p.add_argument("--yrange", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"))
    return p.parse_args()


def run(args: argparse.Namespace, *, verbose: bool = True) -> dict:
    save_outputs = bool(getattr(args, "save_outputs", True))
    write_preview = bool(getattr(args, "write_preview", True)) and save_outputs

    common = prepare_common_inputs(
        args,
        observer_overrides={
            "dsun_cm": args.dsun_cm,
            "lonc_deg": args.lonc_deg,
            "b0sun_deg": args.b0sun_deg,
        },
    )
    model_path = common.model_path
    loader = common.loader

    gxi = gximagecomputing.GXRadioImageComputing()
    model = common.model
    model_dt = common.model_dt
    ebtel_path = common.ebtel_path
    ebtel_c = common.ebtel_c
    ebtel_dt = common.ebtel_dt
    center_source = common.center_source
    xc = common.xc
    yc = common.yc
    dx = common.dx
    dy = common.dy
    nx = common.nx
    ny = common.ny
    fov_x = common.fov_x
    fov_y = common.fov_y

    freqlist = np.arange(5.8, 12.0 + 1e-9, 0.2)[::2].tolist()
    plasma = plasma_defaults()

    result = gxi.synth_model(
        model,
        model_dt,
        ebtel_c,
        ebtel_dt,
        freqlist,
        nx,
        ny,
        xc,
        yc,
        dx,
        dy,
        plasma.tbase,
        plasma.nbase,
        plasma.q0,
        plasma.a,
        plasma.b,
        SHtable=plasma.shtable,
        mode=0,
    )

    out_dir = args.output_dir
    if args.output_name is not None:
        out_name = args.output_name
    else:
        out_name = f"{model_path.name}_py_mw_maps.h5"
    out_path = out_dir / out_name
    obs_time_iso = model_obstime_iso(model)
    sunpy_files = []
    h5_path = None
    preview_path = out_dir / f"{model_path.name}_py_mw_maps_preview.png"
    if save_outputs:
        out_dir.mkdir(parents=True, exist_ok=True)
        if args.output_format in {"h5", "both"}:
            h5_path = save_h5_maps(
                result=result,
                freqlist=freqlist,
                out_h5=out_path,
                model_path=model_path,
                model_format=loader,
                xc=xc,
                yc=yc,
                dx=dx,
                dy=dy,
                obs_time_iso=obs_time_iso,
            )
        if args.output_format in {"fits", "both"}:
            sunpy_files = save_sunpy_maps(
                result=result,
                freqlist=freqlist,
                out_dir=out_dir,
                stem=model_path.stem,
                xc=xc,
                yc=yc,
                dx=dx,
                dy=dy,
                obs_time_iso=obs_time_iso,
            )
        if write_preview:
            if h5_path is not None:
                save_preview(h5_path, preview_path, title=model_path.stem)
            else:
                # Keep preview generation for fits-only mode by writing a transient H5 container.
                preview_h5 = out_dir / f"{model_path.name}_py_mw_maps_preview_tmp.h5"
                save_h5_maps(
                    result=result,
                    freqlist=freqlist,
                    out_h5=preview_h5,
                    model_path=model_path,
                    model_format=loader,
                    xc=xc,
                    yc=yc,
                    dx=dx,
                    dy=dy,
                    obs_time_iso=obs_time_iso,
                )
                save_preview(preview_h5, preview_path, title=model_path.stem)
                preview_h5.unlink(missing_ok=True)

    if verbose:
        print(f"Using library: {gxi.libname}")
        if ebtel_path:
            print(f"Using EBTEL: {ebtel_path}")
        else:
            print("Using EBTEL: none (DEM/heating tables disabled; isothermal/hydrostatic fallback)")
        print(f"Model: {model_path} ({loader})")
        if common.observer_overrides_applied:
            print("Observer overrides: " + ", ".join(f"{k}={v}" for k, v in common.observer_overrides_applied.items()))
        print(f"Center source: {center_source}")
        print(f"Center used: xc={xc:.3f}, yc={yc:.3f} arcsec")
        print(f"FOV={fov_x:.2f}x{fov_y:.2f} arcsec; N={nx}x{ny}; dx={dx:.2f}, dy={dy:.2f} arcsec")
        if save_outputs:
            print("Outputs:")
            if args.output_format in {"h5", "both"} and h5_path is not None:
                print(f"- h5: {h5_path}")
            if args.output_format in {"fits", "both"}:
                print(f"- fits_dir: {out_dir} ({len(sunpy_files)} files)")
                if len(sunpy_files) > 0:
                    print(f"- first_fits: {sunpy_files[0]}")
            if write_preview:
                print(f"- preview_png: {preview_path}")
        else:
            print("Outputs: disabled (save_outputs=False)")

    return {
        "kind": "mw",
        "library_path": str(gxi.libname),
        "model_path": str(model_path),
        "model_format": str(loader),
        "ebtel_path": str(ebtel_path),
        "observer_overrides_applied": dict(common.observer_overrides_applied),
        "center_source": str(center_source),
        "geometry": {
            "xc_arcsec": float(xc),
            "yc_arcsec": float(yc),
            "dx_arcsec": float(dx),
            "dy_arcsec": float(dy),
            "nx": int(nx),
            "ny": int(ny),
            "fov_x_arcsec": float(fov_x),
            "fov_y_arcsec": float(fov_y),
        },
        "obs_time_iso": str(obs_time_iso),
        "freqlist_ghz": [float(f) for f in freqlist],
        "result": result,
        "outputs": {
            "save_outputs": bool(save_outputs),
            "write_preview": bool(write_preview),
            "h5_path": str(h5_path) if h5_path is not None else None,
            "preview_png": (str(preview_path) if write_preview else None),
            "fits_paths": [str(p) for p in sunpy_files],
            "output_dir": str(out_dir),
        },
    }


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
