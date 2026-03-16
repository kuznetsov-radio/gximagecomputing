#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

from gxrender import (
    CoronalPlasmaParameters,
    MapGeometry,
    MWRenderOptions,
    ObserverOverrides,
    render_mw_maps,
)


def _example_default_frequencies() -> list[float]:
    import numpy as np

    return np.arange(5.8, 12.0 + 1e-9, 0.2)[::2].astype(np.float64).tolist()


def _warn_example_default(message: str) -> None:
    warnings.warn(f"SDK MW example default applied: {message}", stacklevel=3)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="SDK example: render MW maps without CLI workflow wrappers.")
    p.add_argument("--model-path", type=Path, required=True, help="Path to CHR model (.sav or .h5)")
    p.add_argument("--model-format", choices=["auto", "sav", "h5"], default="auto")
    p.add_argument("--ebtel-path", type=str, default=None, help="Optional EBTEL table (.sav)")
    p.add_argument("--output-dir", type=Path, default=None, help="Output directory (used only when saving outputs)")
    p.add_argument("--output-name", type=str, default=None, help="Output map filename (H5 or base name for FITS)")
    p.add_argument("--output-format", choices=["h5", "fits", "both"], default="h5")
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--dx", type=float, default=None)
    p.add_argument("--dy", type=float, default=None)
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
    p.add_argument("--pixel-scale-arcsec", type=float, default=None)
    p.add_argument("--frequencies-ghz", nargs="+", type=float, default=None)
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
    p.add_argument("--save-outputs", dest="save_outputs", action="store_true", help="Write outputs (default)")
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
    if args.pixel_scale_arcsec is None and args.dx is None and args.dy is None:
        _warn_example_default("no explicit pixel scale was provided, so dx=dy=2.0 arcsec/pixel was assumed")
        args.pixel_scale_arcsec = 2.0
    if args.frequencies_ghz is None:
        _warn_example_default("no explicit MW frequency list was provided, so 5.8..11.8 GHz was assumed")
        args.frequencies_ghz = _example_default_frequencies()

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

    result = render_mw_maps(
        MWRenderOptions(
            model_path=args.model_path,
            model_format=args.model_format,
            ebtel_path=args.ebtel_path,
            output_dir=args.output_dir,
            output_name=args.output_name,
            output_format=args.output_format,
            freqlist_ghz=args.frequencies_ghz,
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

    print("kind=mw")
    print(f"center_source={result.center_source}")
    print(f"geometry={result.geometry.nx}x{result.geometry.ny} @ {result.geometry.dx_arcsec:.3f}x{result.geometry.dy_arcsec:.3f} arcsec")
    print(f"nfreq={len(result.freqlist_ghz)} first_freqs={result.freqlist_ghz[:3]}")
    print(f"ti shape={tuple(result.ti.shape)}")
    print(f"tv shape={tuple(result.tv.shape)}")
    print(f"save_outputs={result.outputs.save_outputs} write_preview={result.outputs.write_preview}")
    print(f"h5_path={result.outputs.h5_path}")
    print(f"preview_png={result.outputs.preview_png}")
    print(f"n_fits={len(result.outputs.fits_paths)}")


if __name__ == "__main__":
    main()
