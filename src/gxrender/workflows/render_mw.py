from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import gxrender
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.io import fits

from gxrender.io.maps_h5 import save_h5_maps
from gxrender.geometry.observer_geometry import compute_sunpy_wcs_header, observer_summary
from gxrender.workflows._render_common import (
    DEFAULT_OUTDIR,
    model_obstime_iso,
    prepare_common_inputs,
    resolve_mw_frequencies,
    resolve_plasma_parameters,
)


def _example_default_frequencies() -> list[float]:
    return np.arange(5.8, 12.0 + 1e-9, 0.2)[::2].astype(np.float64).tolist()


def _warn_example_default(message: str) -> None:
    warnings.warn(f"Example MW default applied: {message}", stacklevel=3)


def _apply_example_defaults(args: argparse.Namespace) -> argparse.Namespace:
    messages: list[str] = []

    def note(message: str) -> None:
        messages.append(message)

    if args.ebtel_path is None:
        env_ebtel = os.environ.get("GXIMAGECOMPUTING_EBTEL_PATH", "").strip()
        if env_ebtel:
            note(f"no explicit EBTEL path was provided, so GXIMAGECOMPUTING_EBTEL_PATH={env_ebtel!r} was assumed")
            args.ebtel_path = env_ebtel
        else:
            note('no explicit EBTEL path was provided, so "" was used to disable DEM/DDM tables')
            args.ebtel_path = ""

    if args.frequencies_ghz is None:
        note("no explicit MW frequency list was provided, so 5.8..11.8 GHz was assumed")
        args.frequencies_ghz = _example_default_frequencies()

    if args.pixel_scale_arcsec is None and args.dx is None and args.dy is None:
        note("no explicit pixel scale was provided, so dx=dy=2.0 arcsec/pixel was assumed")
        args.pixel_scale_arcsec = 2.0

    plasma_defaults = {
        "tbase": 1e6,
        "nbase": 1e8,
        "q0": 0.0217,
        "a": 0.3,
        "b": 2.7,
    }
    for key, value in plasma_defaults.items():
        if getattr(args, key, None) is None:
            note(f"no explicit {key} was provided, so {value!r} was assumed")
            setattr(args, key, value)

    if args.corona_mode is None:
        note("no explicit corona mode was provided, so mode=0 was assumed")
        args.corona_mode = 0

    if messages:
        _warn_example_default(
            "the following example-layer defaults were applied:\n"
            + "\n".join(f"- {message}" for message in messages)
        )

    return args


def _panel_header(base_wcs_header: fits.Header, **extra_cards: object) -> fits.Header:
    header = base_wcs_header.copy()
    header["NAXIS"] = 2
    for key, value in extra_cards.items():
        header[str(key).upper()] = value
    return header


def _observer_suffix(observer_text: str) -> str:
    compact = str(observer_text).strip()
    if not compact or compact.lower() == "earth":
        return ""
    if compact.lower().startswith("custom"):
        compact = "custom"
    return f" | observer={compact}"


def _preview_header(title: str, freq0: float, date_obs: str, observer_text: str) -> tuple[str, str]:
    top = f"{title} | {freq0:.2f} GHz"
    context = date_obs.strip() if date_obs else ""
    obs_suffix = _observer_suffix(observer_text)
    if obs_suffix:
        context = f"{context}{obs_suffix}" if context else obs_suffix.lstrip(" |")
    return top, context


def save_preview(result: dict, freqlist: list[float], out_png: Path, title: str, base_wcs_header: fits.Header, observer_text: str) -> None:
    ti_cube = np.asarray(result["TI"], dtype=np.float64)
    tv_cube = np.asarray(result["TV"], dtype=np.float64)
    if ti_cube.ndim != 3 or tv_cube.ndim != 3:
        raise ValueError(f"Unexpected rendered-map cube shape: TI={ti_cube.shape}, TV={tv_cube.shape}")
    ti = ti_cube[:, :, 0]
    tv = tv_cube[:, :, 0]
    freq0 = float(freqlist[0]) if len(freqlist) else float("nan")

    m_ti = sunpy.map.Map(ti, _panel_header(base_wcs_header, content="TI", freqghz=freq0))
    m_tv = sunpy.map.Map(tv, _panel_header(base_wcs_header, content="TV", freqghz=freq0))
    date_obs = str(base_wcs_header.get("DATE-OBS", ""))
    header_line, context_line = _preview_header(title, freq0, date_obs, observer_text)

    fig = plt.figure(figsize=(13.6, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[1.0, 0.055, 1.0, 0.055],
        left=0.06,
        right=0.97,
        bottom=0.10,
        top=0.84,
        wspace=0.08,
    )
    ax_ti = fig.add_subplot(gs[0, 0], projection=m_ti)
    ax_tv = fig.add_subplot(gs[0, 2], projection=m_tv)
    im_ti = m_ti.plot(axes=ax_ti, cmap="inferno", interpolation="nearest")
    im_tv = m_tv.plot(axes=ax_tv, cmap="coolwarm", interpolation="nearest")
    fig.suptitle(header_line, fontsize=15, fontweight="semibold", y=0.965)
    if context_line:
        fig.text(0.5, 0.925, context_line, ha="center", va="center", fontsize=10.5, color="#555555")
    ax_ti.set_title("Stokes I", fontsize=12.5, loc="left", pad=8)
    ax_tv.set_title("Stokes V", fontsize=12.5, loc="left", pad=8)
    ax_ti.set_xlabel("Solar X [arcsec]")
    ax_ti.set_ylabel("Solar Y [arcsec]")
    ax_tv.set_xlabel("Solar X [arcsec]")
    ax_tv.set_ylabel("")
    cax_ti = fig.add_subplot(gs[0, 1])
    cax_tv = fig.add_subplot(gs[0, 3])
    cbar0 = fig.colorbar(im_ti, cax=cax_ti, orientation="vertical")
    cbar1 = fig.colorbar(im_tv, cax=cax_tv, orientation="vertical")
    cbar0.set_label("Tb [K]")
    cbar1.set_label("Tb [K]")
    fig.savefig(out_png, dpi=140, facecolor="white")
    plt.close(fig)


def save_sunpy_maps(
    result: dict,
    freqlist: list[float],
    out_dir: Path,
    stem: str,
    base_wcs_header: fits.Header,
) -> list[str]:
    ti = np.asarray(result["TI"], dtype=np.float32)
    tv = np.asarray(result["TV"], dtype=np.float32)
    _ny, _nx, nf = ti.shape
    files: list[str] = []
    for i in range(nf):
        freq = float(freqlist[i])
        m_ti = sunpy.map.Map(ti[:, :, i], _panel_header(base_wcs_header, content="TI", freqghz=freq))
        m_tv = sunpy.map.Map(tv[:, :, i], _panel_header(base_wcs_header, content="TV", freqghz=freq))
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
        help='Optional EBTEL table (.sav). Use "" to explicitly disable DEM/DDM tables.',
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
    p.add_argument("--observer", type=str, default=None, help="Observer name resolved via SunPy (e.g. earth, mars, solo).")
    p.add_argument("--dx", type=float, default=None, help="Pixel scale along X (arcsec/pixel).")
    p.add_argument("--dy", type=float, default=None, help="Pixel scale along Y (arcsec/pixel).")
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
    p.add_argument(
        "--auto-fov",
        action="store_true",
        help="Recompute an inscribing FOV from geometry instead of using saved observer/fov metadata.",
    )
    p.add_argument(
        "--use-saved-fov",
        action="store_true",
        help="Use observer/fov from the input file as the default render window instead of recomputing an inscribing FOV.",
    )
    p.add_argument(
        "--frequencies-ghz",
        "--freqlist-ghz",
        dest="frequencies_ghz",
        type=float,
        nargs="+",
        default=None,
        help="Explicit microwave frequency list in GHz.",
    )
    p.add_argument("--tbase", type=float, default=None, help="Open-field plasma temperature (K).")
    p.add_argument("--nbase", type=float, default=None, help="Open-field base plasma density (cm^-3).")
    p.add_argument("--q0", type=float, default=None, help="Closed-field heating normalization Q0.")
    p.add_argument("--a", type=float, default=None, help="Closed-field heating exponent a.")
    p.add_argument("--b", type=float, default=None, help="Closed-field heating exponent b.")
    p.add_argument(
        "--shtable-path",
        type=Path,
        default=None,
        help="Path to a custom 7x7 selective-heating table (.npy, .npz, .txt, or .csv).",
    )
    p.add_argument(
        "--selective-heating",
        action="store_true",
        help="Enable selective heating. If no SHtable is supplied, the standard 7x7 default table is used.",
    )
    p.add_argument("--corona-mode", type=int, default=None, help="Raw DefineCoronaParams mode bitmask.")
    p.add_argument("--force-isothermal", action="store_true", help="Set DefineCoronaParams /force_isothermal mode bit.")
    p.add_argument("--interpol-b", action="store_true", help="Set DefineCoronaParams interpolB mode bit.")
    p.add_argument("--analytical-nt", action="store_true", help="Set DefineCoronaParams /analyticalNT mode bit.")
    return p.parse_args()


def run(args: argparse.Namespace, *, verbose: bool = True) -> dict:
    save_outputs = bool(getattr(args, "save_outputs", True))
    write_preview = bool(getattr(args, "write_preview", True)) and save_outputs

    common = prepare_common_inputs(args)
    model_path = common.model_path
    loader = common.loader

    gxi = gxrender.GXRadioImageComputing()
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
    observer_geometry = common.observer_geometry
    observer_text = observer_summary(observer_geometry)

    freqlist = resolve_mw_frequencies(args)
    plasma = resolve_plasma_parameters(args)
    warnings.warn(
        "Current Python MW workflow uses projection=0 (parallel off, exact off) for the DLL simbox path. "
        "This assumption is explicit in the workflow today because high-level projection controls are not exposed yet.",
        stacklevel=2,
    )

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
        mode=plasma.mode,
        warn_defaults=False,
    )

    out_dir = args.output_dir
    if args.output_name is not None:
        out_name = args.output_name
    else:
        out_name = f"{model_path.name}_py_mw_maps.h5"
    out_path = out_dir / out_name
    obs_time_iso = model_obstime_iso(model)
    base_wcs_header = compute_sunpy_wcs_header(
        nx=nx,
        ny=ny,
        xc_arcsec=xc,
        yc_arcsec=yc,
        dx_arcsec=dx,
        dy_arcsec=dy,
        obs_time=obs_time_iso,
        observer_geometry=observer_geometry,
        bunit="K",
    )
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
                wcs_header=base_wcs_header,
                observer_name=observer_geometry.observer_name,
                observer_source=observer_geometry.observer_source,
                observer_warnings=observer_geometry.warnings,
                l0_deg=observer_geometry.l0_deg,
                b0_deg=observer_geometry.b0_deg,
                dsun_cm=observer_geometry.dsun_cm,
                rsun_cm=observer_geometry.rsun_cm,
                rsun_arcsec=observer_geometry.rsun_arcsec,
                tbase=plasma.tbase,
                nbase=plasma.nbase,
                q0=plasma.q0,
                a=plasma.a,
                b=plasma.b,
                corona_mode=plasma.mode,
                shtable=plasma.shtable,
            )
        if args.output_format in {"fits", "both"}:
            sunpy_files = save_sunpy_maps(
                result=result,
                freqlist=freqlist,
                out_dir=out_dir,
                stem=model_path.stem,
                base_wcs_header=base_wcs_header,
            )
        if write_preview:
            save_preview(result, freqlist, preview_path, title=model_path.stem, base_wcs_header=base_wcs_header, observer_text=observer_text)

    if verbose:
        print(f"Using library: {gxi.libname}")
        if ebtel_path:
            print(f"Using EBTEL: {ebtel_path}")
        else:
            print("Using EBTEL: none (DEM/heating tables disabled; isothermal/hydrostatic fallback)")
        print(f"Model: {model_path} ({loader})")
        print(f"Observer: {observer_text} ({observer_geometry.observer_source})")
        for warning in observer_geometry.warnings:
            print(f"Observer warning: {warning}")
        print(f"Center source: {center_source}")
        print(f"Center used: xc={xc:.3f}, yc={yc:.3f} arcsec")
        print(f"FOV={fov_x:.2f}x{fov_y:.2f} arcsec; N={nx}x{ny}; dx={dx:.2f}, dy={dy:.2f} arcsec")
        print(
            "MW plasma: "
            f"Tbase={plasma.tbase:.6g} K, nbase={plasma.nbase:.6g} cm^-3, "
            f"Q0={plasma.q0:.6g}, a={plasma.a:.6g}, b={plasma.b:.6g}, mode={plasma.mode}"
        )
        print(
            "SHtable: "
            + (
                "none (selective heating disabled; non-SH DLL entrypoint)"
                if plasma.shtable is None
                else (
                    f"{'custom' if getattr(args, 'shtable_path', None) is not None else ('default' if getattr(args, '_default_shtable_applied', False) else 'provided')} "
                    f"(7x7, min={np.min(plasma.shtable):.6g}, max={np.max(plasma.shtable):.6g})"
                )
            )
        )
        print(
            "MW frequencies [GHz]: "
            f"n={len(freqlist)}, min={min(freqlist):.6g}, max={max(freqlist):.6g}"
        )
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
        "observer": {
            "name": str(observer_geometry.observer_name),
            "source": str(observer_geometry.observer_source),
            "summary": observer_text,
            "warnings": [str(v) for v in observer_geometry.warnings],
            "l0_deg": float(observer_geometry.l0_deg),
            "b0_deg": float(observer_geometry.b0_deg),
            "dsun_cm": float(observer_geometry.dsun_cm),
            "rsun_cm": (float(observer_geometry.rsun_cm) if observer_geometry.rsun_cm is not None else None),
            "rsun_arcsec": (
                float(observer_geometry.rsun_arcsec) if observer_geometry.rsun_arcsec is not None else None
            ),
        },
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
        "plasma": {
            "tbase_k": float(plasma.tbase),
            "nbase_cm3": float(plasma.nbase),
            "q0": float(plasma.q0),
            "a": float(plasma.a),
            "b": float(plasma.b),
            "mode": int(plasma.mode),
            "selective_heating": bool(plasma.selective_heating),
            "shtable": (np.asarray(plasma.shtable, dtype=np.float64).tolist() if plasma.shtable is not None else None),
        },
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
    args = parse_args()
    _apply_example_defaults(args)
    run(args)


if __name__ == "__main__":
    main()
