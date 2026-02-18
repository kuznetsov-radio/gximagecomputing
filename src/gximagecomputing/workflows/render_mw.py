from __future__ import annotations

import argparse
import os
from pathlib import Path

import gximagecomputing
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.time import Time

from gximagecomputing.io.ebtel import resolve_ebtel_path
from gximagecomputing.io.maps_h5 import save_h5_maps
from gximagecomputing.io.model import estimate_hpc_center, infer_center_from_execute, infer_fov_from_execute


# .../gximagecomputing/src/gximagecomputing/workflows/render_mw.py -> repo root at parents[3]
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = (
    Path("C:/Temp/gximagecomputing_validation_groundtruth")
    if os.name == "nt"
    else Path("/tmp/gximagecomputing_validation_groundtruth")
)


def save_preview(result: dict, out_png: Path, title: str, freqlist: list[float], obs_time_iso: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    freq0 = float(freqlist[0]) if len(freqlist) else float("nan")
    im0 = axes[0].imshow(result["TI"][:, :, 0], interpolation=None)
    axes[0].set_title(f"{title} TI @ {freq0:.2f} GHz\\n{obs_time_iso}")
    im1 = axes[1].imshow(result["TV"][:, :, 0], interpolation=None)
    axes[1].set_title(f"{title} TV @ {freq0:.2f} GHz\\n{obs_time_iso}")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    cbar0 = fig.colorbar(im0, ax=axes[0], orientation="vertical", fraction=0.046, pad=0.04)
    cbar1 = fig.colorbar(im1, ax=axes[1], orientation="vertical", fraction=0.046, pad=0.04)
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
        type=Path,
        default=None,
        help="Optional EBTEL table (.sav). If omitted, DEM/heating tables are disabled.",
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


def run(args: argparse.Namespace) -> None:
    os.environ["OMP_NUM_THREADS"] = str(args.omp_threads)

    model_path = args.model_path
    if args.model_format == "auto":
        loader = "h5" if model_path.suffix.lower() in {".h5", ".hdf5"} else "sav"
    else:
        loader = args.model_format

    gxi = gximagecomputing.GXRadioImageComputing()
    ebtel_env = os.environ.get("GXIMAGECOMPUTING_EBTEL_PATH", "").strip()
    if args.ebtel_path is None and not ebtel_env:
        ebtel_path = ""
        ebtel_c, ebtel_dt = gxi.load_ebtel_none()
    else:
        ebtel_path = str(resolve_ebtel_path(args.ebtel_path))
        ebtel_c, ebtel_dt = gxi.load_ebtel(ebtel_path)

    if loader == "h5":
        model, model_dt = gxi.load_model_hdf(str(model_path))
    else:
        model, model_dt = gxi.load_model_sav(str(model_path))

    center_exec = infer_center_from_execute(loader_name=loader, model_path=model_path)
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
        if args.nx is not None:
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
        if args.ny is not None:
            print("Note: --yrange provided, so --ny is ignored and recomputed from --dy.")
        yc = 0.5 * (ymin + ymax)
        ny = max(16, int(np.ceil((ymax - ymin) / dy)))
    else:
        yc = float(args.yc) if args.yc is not None else float(yc_auto)
        ny = max(16, int(args.ny)) if args.ny is not None else max(16, int(np.ceil(model_h_arcsec / dy)))

    fov_x = nx * dx
    fov_y = ny * dy

    freqlist = np.arange(5.8, 12.0 + 1e-9, 0.2)[::2].tolist()
    # Coronal plasma parameters (aligned with RenderExampleMW.pro defaults).
    tbase = 1e6
    nbase = 1e8
    q0 = 0.0217
    a = 0.3
    b = 2.7
    w = np.array([1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4], dtype=np.float64)
    shtable = np.outer(w, w)
    shtable[6, 6] = 0.1

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
        tbase,
        nbase,
        q0,
        a,
        b,
        SHtable=shtable,
        mode=0,
    )

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.output_name is not None:
        out_name = args.output_name
    else:
        out_name = f"{model_path.name}_py_mw_maps.h5"
    out_path = out_dir / out_name
    obs_time_iso = Time(float(model["obstime"][0]) + 283996800.0, format="unix").isot
    sunpy_files = []
    h5_path = None
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
    preview_path = out_dir / f"{model_path.name}_py_mw_maps_preview.png"
    save_preview(result, preview_path, title=model_path.stem, freqlist=freqlist, obs_time_iso=obs_time_iso)

    print(f"Using library: {gxi.libname}")
    if ebtel_path:
        print(f"Using EBTEL: {ebtel_path}")
    else:
        print("Using EBTEL: none (DEM/heating tables disabled; isothermal/hydrostatic fallback)")
    print(f"Model: {model_path} ({loader})")
    print(f"Center source: {center_source}")
    print(f"Center used: xc={xc:.3f}, yc={yc:.3f} arcsec")
    print(f"FOV={fov_x:.2f}x{fov_y:.2f} arcsec; N={nx}x{ny}; dx={dx:.2f}, dy={dy:.2f} arcsec")
    print("Outputs:")
    if args.output_format in {"h5", "both"} and h5_path is not None:
        print(f"- h5: {h5_path}")
    if args.output_format in {"fits", "both"}:
        print(f"- fits_dir: {out_dir} ({len(sunpy_files)} files)")
        if len(sunpy_files) > 0:
            print(f"- first_fits: {sunpy_files[0]}")
    print(f"- preview_png: {preview_path}")


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
