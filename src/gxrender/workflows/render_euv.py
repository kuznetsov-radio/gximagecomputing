from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.io import fits

from gxrender.euv import GXEUVImageComputing, build_default_aia_euv_response, load_euv_response_sav
from gxrender.geometry.observer_geometry import compute_sunpy_wcs_header, observer_summary
from gxrender.io.maps_h5 import save_h5_euv_maps
from gxrender.policy.contracts import EUVResponseRequest
from gxrender.policy.euv_response_policy import apply_default_response_selection, resolve_euv_response
from gxrender.utils.render_map_view import _euv_sunpy_cmap_name
from gxrender.utils.test_data import test_data_setup_hint, try_find_response_file
from gxrender.workflows._render_common import (
    DEFAULT_OUTDIR,
    model_obstime_iso,
    prepare_common_inputs,
    resolve_plasma_parameters,
)


def _warn_example_default(message: str) -> None:
    warnings.warn(f"Example EUV default applied: {message}", stacklevel=3)


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


def _resolve_default_response_sav(instrument: str) -> Path | None:
    env_path = os.environ.get("GXIMAGECOMPUTING_EUV_RESPONSE_SAV", "").strip()
    if env_path:
        p = Path(env_path).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"GXIMAGECOMPUTING_EUV_RESPONSE_SAV points to a missing file: {p}")

    return try_find_response_file(str(instrument).strip().lower())


def _apply_default_response_selection(args: argparse.Namespace, *, observer_geometry) -> argparse.Namespace:
    if getattr(args, "response", None) is not None or args.response_sav is not None or args.instrument is not None:
        return args

    observer_name = str(getattr(observer_geometry, "observer_name", "") or "").strip().lower()
    observer_source = str(getattr(observer_geometry, "observer_source", "") or "").strip().lower()

    if observer_name == "stereo-a":
        args.instrument = "STEREO-A"
        return args
    if observer_name == "stereo-b":
        args.instrument = "STEREO-B"
        return args
    if observer_name == "solar orbiter":
        args.instrument = "SOLO-FSI"
        return args
    if observer_name == "earth" and observer_source != "default_earth":
        args.instrument = "AIA"
        return args
    if observer_source == "default_earth":
        _warn_example_default(
            "no explicit instrument name or observer metadata was provided, so AIA was assumed"
        )
        args.instrument = "AIA"
        return args

    raise ValueError(
        "No explicit EUV instrument or response was provided. "
        f"Observer metadata resolves to {observer_name or '<unknown>'!r}, but no default EUV instrument mapping is defined. "
        "Provide --instrument or --response-sav explicitly."
    )


def _resolve_response_inputs(args: argparse.Namespace, *, obs_time_iso: str):
    if getattr(args, "response", None) is not None:
        missing = [
            name
            for name in ("response_dt", "response_meta")
            if getattr(args, name, None) is None
        ]
        if missing:
            raise ValueError(
                "Prebuilt EUV response inputs are incomplete. "
                f"Missing: {', '.join(missing)}."
            )
        return args.response, args.response_dt, args.response_meta

    if args.response_sav is not None:
        return load_euv_response_sav(str(args.response_sav))

    provider_error = None
    if str(args.instrument).strip().lower() == "aia":
        try:
            return build_default_aia_euv_response(
                obstime=obs_time_iso,
                channels=None if args.channels is None else [str(channel) for channel in args.channels],
                correction_state="evenorm_chiantifix",
            )
        except (ImportError, FileNotFoundError) as exc:
            provider_error = exc

    auto_response_sav = _resolve_default_response_sav(args.instrument)
    if auto_response_sav is not None:
        return load_euv_response_sav(str(auto_response_sav))

    hint = test_data_setup_hint(f"EUV response file for instrument {str(args.instrument).strip().lower()!r}")
    if provider_error is not None:
        raise FileNotFoundError(
            "No explicit EUV response was provided. "
            "The Python-native AIA response provider was unavailable, and no default response fixture could be found. "
            f"Provider error: {provider_error}. {hint} "
            "You may also set GXIMAGECOMPUTING_EUV_RESPONSE_SAV to an explicit SAV file."
        ) from provider_error
    raise FileNotFoundError(
        "No explicit EUV response was provided, and no default response fixture could be found. "
        + hint
        + " You may also set GXIMAGECOMPUTING_EUV_RESPONSE_SAV to an explicit SAV file."
    )


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


def _preview_header(title: str, instrument: str, channel: str, date_obs: str, observer_text: str) -> tuple[str, str]:
    inst = str(instrument).strip() or "EUV"
    top = f"{title} | {inst} {channel}"
    context = date_obs.strip() if date_obs else ""
    obs_suffix = _observer_suffix(observer_text)
    if obs_suffix:
        context = f"{context}{obs_suffix}" if context else obs_suffix.lstrip(" |")
    return top, context


def _apply_log_scaling(flux: np.ndarray) -> np.ndarray:
    """Apply log10 scaling to enhance dynamic range, handling invalid values."""
    flux = np.asarray(flux, dtype=np.float64)
    # Mask invalid values
    valid = (np.isfinite(flux)) & (flux > 0)
    if not valid.any():
        return np.ones_like(flux)  # Return ones if all values are invalid
    min_valid = np.nanmin(flux[valid])
    # Replace invalid/negative with min valid value, then apply log10
    flux_safe = np.where(valid, flux, min_valid)
    return np.log10(flux_safe)


def _preview_euv(
    flux_cor: np.ndarray,
    flux_tr: np.ndarray,
    instrument: str,
    channel: str,
    out_png: Path | None,
    title: str,
    base_wcs_header: fits.Header,
    observer_text: str,
    show: bool = False,
    log_scale: bool = False,
):
    if log_scale:
        flux_cor = _apply_log_scaling(flux_cor)
        flux_tr = _apply_log_scaling(flux_tr)
    m_cor = sunpy.map.Map(flux_cor, _panel_header(base_wcs_header, content="EUV_CORONA", channel=str(channel)))
    m_tr = sunpy.map.Map(flux_tr, _panel_header(base_wcs_header, content="EUV_TR", channel=str(channel)))
    m_sum = sunpy.map.Map((flux_cor + flux_tr), _panel_header(base_wcs_header, content="EUV_SUM", channel=str(channel)))
    date_obs = str(base_wcs_header.get("DATE-OBS", ""))
    header_line, context_line = _preview_header(title, str(instrument), str(channel), date_obs, observer_text)
    channel_cmap = _euv_sunpy_cmap_name(str(instrument), str(channel)) or "magma"

    fig = plt.figure(figsize=(16.2, 5.1), constrained_layout=True)
    gs = fig.add_gridspec(
        1,
        6,
        width_ratios=[1.0, 0.05, 1.0, 0.05, 1.0, 0.05],
        left=0.05,
        right=0.975,
        bottom=0.10,
        top=0.84,
        wspace=0.10,
    )
    ax1 = fig.add_subplot(gs[0, 0], projection=m_cor)
    ax2 = fig.add_subplot(gs[0, 2], projection=m_tr)
    ax3 = fig.add_subplot(gs[0, 4], projection=m_sum)
    im1 = m_cor.plot(axes=ax1, cmap=channel_cmap, interpolation="nearest")
    im2 = m_tr.plot(axes=ax2, cmap=channel_cmap, interpolation="nearest")
    im3 = m_sum.plot(axes=ax3, cmap=channel_cmap, interpolation="nearest")
    fig.suptitle(header_line, fontsize=15, fontweight="semibold", y=0.965)
    if context_line:
        fig.text(0.5, 0.925, context_line, ha="center", va="center", fontsize=10.5, color="#555555")
    ax1.set_title("Corona", fontsize=12.5, loc="left", pad=8)
    ax2.set_title("Transition Region", fontsize=12.5, loc="left", pad=8)
    ax3.set_title("Total", fontsize=12.5, loc="left", pad=8)
    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("Solar X [arcsec]")
        ax.set_ylabel("Solar Y [arcsec]")
    ax2.set_ylabel("")
    ax3.set_ylabel("")
    fig.colorbar(im1, cax=fig.add_subplot(gs[0, 1]), orientation="vertical").set_label("DN s^-1 pix^-1")
    fig.colorbar(im2, cax=fig.add_subplot(gs[0, 3]), orientation="vertical").set_label("DN s^-1 pix^-1")
    fig.colorbar(im3, cax=fig.add_subplot(gs[0, 5]), orientation="vertical").set_label("DN s^-1 pix^-1")
    if out_png is not None:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=140, facecolor="white")
    if show:
        plt.show()
    plt.close(fig)


def _resolve_preview_channel_index(channels: list[str], configured_index: int | None) -> int:
    if not channels:
        return 0
    index_value = int(configured_index if configured_index is not None else 0)
    return int(np.clip(index_value, 0, len(channels) - 1))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render EUV maps from a single H5/SAV CHR model.")
    p.add_argument("--model-path", type=Path, required=True, help="Path to input CHR model (.h5 or .sav).")
    p.add_argument("--model-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument(
        "--ebtel-path",
        type=str,
        default=None,
        help='Optional EBTEL table (.sav). Use "" to explicitly disable DEM/DDM tables.',
    )
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--output-name", type=str, default=None, help="Default: <model-basename>_py_euv_maps.h5")
    p.add_argument("--channels", nargs="*", default=None)
    p.add_argument("--instrument", type=str, default=None)
    p.add_argument(
        "--response-sav",
        type=Path,
        default=None,
        help="Path to an IDL gxresponse SAV (for example resp_aia_20251126T153431.sav).",
    )
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--dsun-cm", type=float, default=None, help="Override model.DSun before rendering (cm).")
    p.add_argument("--lonc-deg", type=float, default=None, help="Override model.lonC before rendering (deg).")
    p.add_argument("--b0sun-deg", type=float, default=None, help="Override model.b0Sun before rendering (deg).")
    p.add_argument("--observer", type=str, default=None, help="Observer name resolved via SunPy (e.g. earth, mars, solo).")
    p.add_argument("--dx", type=float, default=None)
    p.add_argument("--dy", type=float, default=None)
    p.add_argument("--pixel-scale-arcsec", type=float, default=None)
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
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
    p.add_argument(
        "--preview-channel-index",
        type=int,
        default=0,
        help="Channel index used for the workflow-generated preview PNG (default: 0).",
    )
    p.add_argument(
        "--preview-png",
        type=Path,
        default=None,
        help="Optional explicit output path for preview PNG. Default: <model>_py_euv_maps_preview_<channel>.png",
    )
    p.add_argument(
        "--log-scale",
        action="store_true",
        help="Apply log10 scaling to preview PNG for enhanced dynamic range visualization.",
    )
    return p.parse_args()

def run(args: argparse.Namespace, *, verbose: bool = True) -> dict:
    save_outputs = bool(getattr(args, "save_outputs", True))
    write_preview = bool(getattr(args, "write_preview", True)) and save_outputs

    # EUV parity is better when the default center is derived from model WCS
    # metadata (IDL-like path) rather than the execute-string shortcut.
    common = prepare_common_inputs(args, prefer_execute_center=False)
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
    observer_geometry = common.observer_geometry
    observer_text = observer_summary(observer_geometry)

    apply_default_response_selection(args, observer_geometry=observer_geometry)

    resolved_response = resolve_euv_response(EUVResponseRequest(args=args, obs_time_iso=model_obstime_iso(model)))
    response = resolved_response.response
    response_dt = resolved_response.response_dt
    response_meta = resolved_response.response_meta

    plasma = resolve_plasma_parameters(args)
    warnings.warn(
        "Current Python EUV workflow uses projection flags off (parallel=False, exact=False, nthreads=0) for the DLL simbox path. "
        "This assumption is explicit in the workflow today because high-level projection controls are not exposed yet.",
        stacklevel=2,
    )

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
        mode=plasma.mode,
        shtable=plasma.shtable,
        warn_defaults=False,
    )

    out_dir = args.output_dir
    out_name = args.output_name if args.output_name is not None else f"{model_path.name}_py_euv_maps.h5"
    out_h5 = out_dir / out_name
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
        bunit="DN s^-1 pix^-1",
    )
    preview_idx = _resolve_preview_channel_index([str(c) for c in response_meta.channels], getattr(args, "preview_channel_index", 0))
    preview_channel = str(response_meta.channels[preview_idx]) if response_meta.channels else ""
    preview_path = (
        Path(args.preview_png)
        if getattr(args, "preview_png", None) is not None
        else out_dir / f"{model_path.name}_py_euv_maps_preview_{preview_channel}.png"
    )
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
        if write_preview:
            _preview_euv(
                flux_cor=out["flux_corona"][:, :, preview_idx],
                flux_tr=out["flux_tr"][:, :, preview_idx],
                instrument=response_meta.instrument,
                channel=preview_channel,
                out_png=preview_path,
                title=model_path.stem,
                base_wcs_header=base_wcs_header,
                observer_text=observer_text,
                log_scale=bool(getattr(args, "log_scale", False)),
            )

    if verbose:
        print(f"Using library: {gx.libname}")
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
            "Coronal plasma: "
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
        print(f"Response: {response_meta.instrument} channels={','.join(response_meta.channels)}")
        if response_meta.source:
            print(f"Response source: {response_meta.source} ({response_meta.mode})")
        else:
            print(f"Response source: prebuilt response object ({response_meta.mode})")
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
    args = parse_args()
    _apply_example_defaults(args)
    run(args)


if __name__ == "__main__":
    main()
