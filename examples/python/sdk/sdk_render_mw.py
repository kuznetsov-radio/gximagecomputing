#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from gximagecomputing import MapGeometry, MWRenderOptions, ObserverOverrides, render_mw_maps


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
    result = render_mw_maps(
        MWRenderOptions(
            model_path=args.model_path,
            model_format=args.model_format,
            ebtel_path=args.ebtel_path,
            output_dir=args.output_dir,
            output_name=args.output_name,
            output_format=args.output_format,
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
