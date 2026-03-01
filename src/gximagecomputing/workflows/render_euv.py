from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map

from gximagecomputing.euv import GXEUVImageComputing, build_default_euv_response, load_euv_response_sav
from gximagecomputing.io.maps_h5 import save_h5_euv_maps
from gximagecomputing.workflows._render_common import DEFAULT_OUTDIR, model_obstime_iso, plasma_defaults, prepare_common_inputs


def _resolve_default_response_sav(instrument: str) -> Path | None:
    env_path = os.environ.get("GXIMAGECOMPUTING_EUV_RESPONSE_SAV", "").strip()
    if env_path:
        p = Path(env_path).expanduser()
        if p.exists():
            return p

    ssw = os.environ.get("SSW", "").strip()
    if not ssw:
        return None

    base = Path(ssw).expanduser() / "packages" / "gx_simulator" / "euv"
    inst = str(instrument).strip().lower()
    name_by_inst = {
        "aia": "aia_response.sav",
        "aia2": "aia_response.sav",
        "trace": "trace_response.sav",
        "sxt": "sxt_response.sav",
    }
    fname = name_by_inst.get(inst)
    if fname is None:
        return None
    for p in (base / fname, base / "response" / fname, base / "responses" / fname):
        if p.exists():
            return p
    return None


def _preview_euv(flux_cor: np.ndarray, flux_tr: np.ndarray, channel: str, out_png: Path, title: str, obs_time_iso: str, xc: float, yc: float, dx: float, dy: float):
    ny, nx = flux_cor.shape[0], flux_cor.shape[1]
    meta = {
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
        "bunit": "DN s^-1 pix^-1",
    }
    m_cor = sunpy.map.Map(flux_cor, dict(meta, content="EUV_CORONA", channel=str(channel)))
    m_tr = sunpy.map.Map(flux_tr, dict(meta, content="EUV_TR", channel=str(channel)))
    m_sum = sunpy.map.Map((flux_cor + flux_tr), dict(meta, content="EUV_SUM", channel=str(channel)))

    fig = plt.figure(figsize=(12, 3.8))
    ax1 = fig.add_subplot(1, 3, 1, projection=m_cor)
    ax2 = fig.add_subplot(1, 3, 2, projection=m_tr)
    ax3 = fig.add_subplot(1, 3, 3, projection=m_sum)
    im1 = m_cor.plot(axes=ax1, cmap="magma", interpolation="nearest")
    im2 = m_tr.plot(axes=ax2, cmap="magma", interpolation="nearest")
    im3 = m_sum.plot(axes=ax3, cmap="magma", interpolation="nearest")
    ax1.set_title(f"{title} COR {channel}\\n{obs_time_iso}")
    ax2.set_title(f"{title} TR {channel}\\n{obs_time_iso}")
    ax3.set_title(f"{title} COR+TR {channel}\\n{obs_time_iso}")
    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("Solar X [arcsec]")
        ax.set_ylabel("Solar Y [arcsec]")
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04).set_label("DN s^-1 pix^-1")
    fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04).set_label("DN s^-1 pix^-1")
    fig.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04).set_label("DN s^-1 pix^-1")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render EUV maps from a single H5/SAV CHR model.")
    p.add_argument("--model-path", type=Path, required=True, help="Path to input CHR model (.h5 or .sav).")
    p.add_argument("--model-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument(
        "--ebtel-path",
        type=str,
        default=None,
        help='Optional EBTEL table (.sav). Use "" to explicitly disable and force isothermal mode.',
    )
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--output-name", type=str, default=None, help="Default: <model-basename>_py_euv_maps.h5")
    p.add_argument("--channels", nargs="*", default=["94", "131", "171", "193", "211", "304", "335"])
    p.add_argument("--instrument", type=str, default="AIA")
    p.add_argument("--response-sav", type=Path, default=None, help="Path to IDL gxresponse SAV (e.g., aia_response.sav)")
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--dsun-cm", type=float, default=None, help="Override model.DSun before rendering (cm).")
    p.add_argument("--lonc-deg", type=float, default=None, help="Override model.lonC before rendering (deg).")
    p.add_argument("--b0sun-deg", type=float, default=None, help="Override model.b0Sun before rendering (deg).")
    p.add_argument("--dx", type=float, default=None)
    p.add_argument("--dy", type=float, default=None)
    p.add_argument("--pixel-scale-arcsec", type=float, default=None)
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
    p.add_argument("--xrange", type=float, nargs=2, default=None, metavar=("XMIN", "XMAX"))
    p.add_argument("--yrange", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"))
    return p.parse_args()


def run(args: argparse.Namespace, *, verbose: bool = True) -> dict:
    save_outputs = bool(getattr(args, "save_outputs", True))
    write_preview = bool(getattr(args, "write_preview", True)) and save_outputs

    # EUV parity is better when the default center is derived from model WCS
    # metadata (IDL-like path) rather than the execute-string shortcut.
    common = prepare_common_inputs(
        args,
        prefer_execute_center=False,
        observer_overrides={
            "dsun_cm": args.dsun_cm,
            "lonc_deg": args.lonc_deg,
            "b0sun_deg": args.b0sun_deg,
        },
    )
    model_path = common.model_path
    loader = common.loader
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

    if args.response_sav is not None:
        response, response_dt, response_meta = load_euv_response_sav(str(args.response_sav))
    else:
        auto_response_sav = _resolve_default_response_sav(args.instrument)
        if auto_response_sav is not None:
            response, response_dt, response_meta = load_euv_response_sav(str(auto_response_sav))
        else:
            response, response_dt, response_meta = build_default_euv_response(
                instrument=args.instrument, channels=[str(c) for c in args.channels]
            )

    plasma = plasma_defaults()

    gx = GXEUVImageComputing()
    out = gx.synth_euv(
        model=model,
        model_dt=model_dt,
        ebtel=ebtel_c,
        ebtel_dt=ebtel_dt,
        response=response,
        response_dt=response_dt,
        box_nx=nx,
        box_ny=ny,
        box_xc=xc,
        box_yc=yc,
        box_dx=dx,
        box_dy=dy,
        tbase=plasma.tbase,
        nbase=plasma.nbase,
        q0=plasma.q0,
        a=plasma.a,
        b=plasma.b,
        mode=0,
        shtable=plasma.shtable,
    )

    out_dir = args.output_dir
    out_name = args.output_name if args.output_name is not None else f"{model_path.name}_py_euv_maps.h5"
    out_h5 = out_dir / out_name
    obs_time_iso = model_obstime_iso(model)
    preview_path = out_dir / f"{model_path.name}_py_euv_maps_preview.png"
    h5_path = None
    if save_outputs:
        out_dir.mkdir(parents=True, exist_ok=True)
        h5_path = save_h5_euv_maps(
            flux_corona=out["flux_corona"],
            flux_tr=out["flux_tr"],
            channels=response_meta.channels,
            out_h5=out_h5,
            model_path=model_path,
            model_format=loader,
            xc=xc,
            yc=yc,
            dx=dx,
            dy=dy,
            obs_time_iso=obs_time_iso,
            instrument=response_meta.instrument,
        )
        if write_preview:
            _preview_euv(
                flux_cor=out["flux_corona"][:, :, 0],
                flux_tr=out["flux_tr"][:, :, 0],
                channel=response_meta.channels[0],
                out_png=preview_path,
                title=model_path.stem,
                obs_time_iso=obs_time_iso,
                xc=xc,
                yc=yc,
                dx=dx,
                dy=dy,
            )

    if verbose:
        print(f"Using library: {gx.libname}")
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
        print(f"Response: {response_meta.instrument} channels={','.join(response_meta.channels)}")
        if response_meta.source:
            print(f"Response source: {response_meta.source} ({response_meta.mode})")
        else:
            print("Response source: synthetic fallback (not IDL-equivalent); pass --response-sav or set GXIMAGECOMPUTING_EUV_RESPONSE_SAV")
        if save_outputs:
            print("Outputs:")
            print(f"- h5: {h5_path}")
            if write_preview:
                print(f"- preview_png: {preview_path}")
        else:
            print("Outputs: disabled (save_outputs=False)")

    return {
        "kind": "euv",
        "library_path": str(gx.libname),
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
        "response": {
            "instrument": str(response_meta.instrument),
            "channels": [str(c) for c in response_meta.channels],
            "source": str(response_meta.source),
            "mode": str(response_meta.mode),
        },
        "result": out,
        "outputs": {
            "save_outputs": bool(save_outputs),
            "write_preview": bool(write_preview),
            "h5_path": (str(h5_path) if h5_path is not None else None),
            "preview_png": (str(preview_path) if write_preview else None),
            "output_dir": str(out_dir),
        },
    }


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
