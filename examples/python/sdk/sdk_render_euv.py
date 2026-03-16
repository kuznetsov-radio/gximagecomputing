#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

from gxrender import (
    CoronalPlasmaParameters,
    EUVRenderOptions,
    MapGeometry,
    ObserverOverrides,
    render_euv_maps,
)
from gxrender.utils.test_data import test_data_setup_hint, try_find_response_file


def _warn_example_default(message: str) -> None:
    warnings.warn(f"SDK EUV example default applied: {message}", stacklevel=3)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SDK example: render EUV maps without CLI workflow wrappers.")
    p.add_argument("--model-path", type=Path, required=True, help="Path to CHR model (.sav or .h5)")
    p.add_argument("--model-format", choices=["auto", "sav", "h5"], default="auto")
    p.add_argument("--ebtel-path", type=str, default=None, help="Optional EBTEL table (.sav)")
    p.add_argument("--response-sav", type=Path, default=None, help="Optional EUV response SAV file")
    p.add_argument("--output-dir", type=Path, default=None, help="Output directory (used only when saving outputs)")
    p.add_argument("--output-name", type=str, default=None, help="Output H5 filename")
    p.add_argument("--instrument", type=str, default=None)
    p.add_argument("--channels", nargs="+", default=None, help="Channel list (e.g. 94 131 171 193 211 304 335)")
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--dx", type=float, default=None)
    p.add_argument("--dy", type=float, default=None)
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
    p.add_argument("--pixel-scale-arcsec", type=float, default=None)
    p.add_argument("--tbase", type=float, default=None)
    p.add_argument("--nbase", type=float, default=None)
    p.add_argument("--q0", type=float, default=None)
    p.add_argument("--a", type=float, default=None)
    p.add_argument("--b", type=float, default=None)
    p.add_argument("--corona-mode", type=int, default=None)
    p.add_argument("--shtable-path", type=Path, default=None)
    p.add_argument("--selective-heating", action="store_true")
    p.add_argument("--dsun-cm", type=float, default=None)
    p.add_argument("--lonc-deg", type=float, default=None)
    p.add_argument("--b0sun-deg", type=float, default=None)
    p.add_argument("--save-outputs", dest="save_outputs", action="store_true", help="Write H5 outputs (default)")
    p.add_argument("--no-save-outputs", dest="save_outputs", action="store_false", help="Keep render in memory only")
    p.set_defaults(save_outputs=True)
    p.add_argument("--write-preview", dest="write_preview", action="store_true", help="Write preview PNG (default)")
    p.add_argument("--no-write-preview", dest="write_preview", action="store_false", help="Disable preview PNG")
    p.set_defaults(write_preview=True)
    p.add_argument("--verbose", action="store_true", help="Enable workflow progress logging")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.ebtel_path is None:
        env_ebtel = os.environ.get("GXIMAGECOMPUTING_EBTEL_PATH", "").strip()
        if env_ebtel:
            _warn_example_default(f"no explicit EBTEL path was provided, so GXIMAGECOMPUTING_EBTEL_PATH={env_ebtel!r} was assumed")
            args.ebtel_path = env_ebtel
        else:
            _warn_example_default('no explicit EBTEL path was provided, so "" was used to disable DEM/DDM tables')
            args.ebtel_path = ""
    if args.instrument is None:
        _warn_example_default("no explicit instrument name was provided, so AIA was assumed")
        args.instrument = "AIA"
    if args.channels is None:
        _warn_example_default("no explicit EUV channel list was provided, so the standard AIA channels were assumed")
        args.channels = ["94", "131", "171", "193", "211", "304", "335"]
    if args.pixel_scale_arcsec is None and args.dx is None and args.dy is None:
        _warn_example_default("no explicit pixel scale was provided, so dx=dy=2.0 arcsec/pixel was assumed")
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
            _warn_example_default(f"no explicit {key} was provided, so {value!r} was assumed")
            setattr(args, key, value)
    if args.corona_mode is None:
        _warn_example_default("no explicit corona mode was provided, so mode=0 was assumed")
        args.corona_mode = 0

    shtable = None
    if args.shtable_path is not None:
        import numpy as np

        shtable = np.load(args.shtable_path)

    if args.response_sav is None:
        env_response = os.environ.get("GXIMAGECOMPUTING_EUV_RESPONSE_SAV", "").strip()
        if env_response:
            response_path = Path(env_response).expanduser()
            if not response_path.exists():
                raise FileNotFoundError(
                    f"GXIMAGECOMPUTING_EUV_RESPONSE_SAV points to a missing file: {response_path}"
                )
            _warn_example_default(
                f"no explicit EUV response SAV was provided, so GXIMAGECOMPUTING_EUV_RESPONSE_SAV={env_response!r} was assumed"
            )
            args.response_sav = response_path
        else:
            auto_response_sav = try_find_response_file(args.instrument)
            if auto_response_sav is None:
                raise FileNotFoundError(
                    "No explicit EUV response SAV was provided, and no default response fixture could be found. "
                    + test_data_setup_hint(f"EUV response file for instrument {str(args.instrument).strip().lower()!r}")
                    + " You may also set GXIMAGECOMPUTING_EUV_RESPONSE_SAV to an explicit SAV file."
                )
            _warn_example_default(
                f"no explicit EUV response SAV was provided, so {auto_response_sav} was assumed"
            )
            args.response_sav = auto_response_sav

    response = None
    response_dt = None
    response_meta = None

    result = render_euv_maps(
        EUVRenderOptions(
            model_path=args.model_path,
            model_format=args.model_format,
            ebtel_path=args.ebtel_path,
            output_dir=args.output_dir,
            output_name=args.output_name,
            channels=args.channels,
            instrument=args.instrument,
            response_sav=args.response_sav,
            response=response,
            response_dt=response_dt,
            response_meta=response_meta,
            plasma=CoronalPlasmaParameters(
                tbase=args.tbase,
                nbase=args.nbase,
                q0=args.q0,
                a=args.a,
                b=args.b,
                mode=args.corona_mode,
                selective_heating=args.selective_heating,
                shtable=shtable,
            ),
            omp_threads=args.omp_threads,
            geometry=MapGeometry(
                xc=args.xc,
                yc=args.yc,
                dx=args.dx,
                dy=args.dy,
                nx=args.nx,
                ny=args.ny,
                pixel_scale_arcsec=args.pixel_scale_arcsec,
            ),
            observer=ObserverOverrides(
                dsun_cm=args.dsun_cm,
                lonc_deg=args.lonc_deg,
                b0sun_deg=args.b0sun_deg,
            ),
            save_outputs=args.save_outputs,
            write_preview=args.write_preview,
            verbose=args.verbose,
        )
    )

    print("kind=euv")
    print(f"center_source={result.center_source}")
    print(f"geometry={result.geometry.nx}x{result.geometry.ny} @ {result.geometry.dx_arcsec:.3f}x{result.geometry.dy_arcsec:.3f} arcsec")
    print(f"response={result.response.instrument} channels={','.join(result.response.channels)}")
    print(f"flux_corona shape={tuple(result.flux_corona.shape)}")
    print(f"flux_tr shape={tuple(result.flux_tr.shape)}")
    print(f"save_outputs={result.outputs.save_outputs} write_preview={result.outputs.write_preview}")
    print(f"h5_path={result.outputs.h5_path}")
    print(f"preview_png={result.outputs.preview_png}")


if __name__ == "__main__":
    main()
