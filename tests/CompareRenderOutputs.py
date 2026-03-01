#!/usr/bin/env python3

from __future__ import annotations

# Runtime notes:
# - If running from repository source (without pip install), use:
#     PYTHONPATH=src python tests/CompareRenderOutputs.py ...
# - If SunPy/Matplotlib default config/cache folders are not writable, use:
#     SUNPY_CONFIGDIR=/tmp/sunpy_cfg MPLCONFIGDIR=/tmp/mpl_cfg python ...

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
RENDER_SCRIPT = REPO_ROOT / "examples" / "python" / "cli" / "RenderExampleMW.py"
DEFAULT_TMP_DIR = Path(tempfile.gettempdir())


def _percent_diff_vs_model2(model1: np.ndarray, model2: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    out = np.full(model1.shape, np.nan, dtype=np.float64)
    denom = np.abs(model2)
    m = denom > eps
    out[m] = 100.0 * (model1[m] - model2[m]) / denom[m]
    return out


def _panel_4x4(cube: np.ndarray, freqs: np.ndarray, out_png: Path, title_prefix: str) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(14, 12))
    for i, ax in enumerate(axes.flat):
        if i >= cube.shape[-1]:
            ax.axis("off")
            continue
        im = ax.imshow(cube[:, :, i], interpolation=None)
        ax.set_title(f"{title_prefix} {freqs[i]:.1f} GHz", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _stats(arr: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.nanmean(arr)),
        "std": float(np.nanstd(arr)),
        "min": float(np.nanmin(arr)),
        "max": float(np.nanmax(arr)),
    }


def _run_render(
    model_path: Path,
    model_format: str,
    output_dir: Path,
    output_name: str,
    ebtel_path: Path | None,
    xc: float | None,
    yc: float | None,
    nx: int | None,
    ny: int | None,
    pixel_scale_arcsec: float,
    omp_threads: int,
) -> Path:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")

    cmd = [
        sys.executable,
        str(RENDER_SCRIPT),
        "--model-path",
        str(model_path),
        "--model-format",
        model_format,
        "--output-dir",
        str(output_dir),
        "--output-name",
        output_name,
        "--output-format",
        "h5",
        "--pixel-scale-arcsec",
        str(pixel_scale_arcsec),
        "--omp-threads",
        str(omp_threads),
    ]
    if ebtel_path is not None:
        cmd.extend(["--ebtel-path", str(ebtel_path)])
    if xc is not None:
        cmd.extend(["--xc", str(xc)])
    if yc is not None:
        cmd.extend(["--yc", str(yc)])
    if nx is not None:
        cmd.extend(["--nx", str(nx)])
    if ny is not None:
        cmd.extend(["--ny", str(ny)])

    subprocess.run(cmd, check=True, env=env)
    return output_dir / output_name


def _load_render_h5(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # (nx, ny, nf, 2)
        freqs = np.asarray(f["maps"]["freqlist_ghz"], dtype=np.float64)
        xc = float(f["metadata"]["xc_arcsec"][()])
        yc = float(f["metadata"]["yc_arcsec"][()])
        nx = int(f["metadata"]["nx"][()])
        ny = int(f["metadata"]["ny"][()])
    if cube.ndim != 4 or cube.shape[-1] != 2:
        raise ValueError(f"Unexpected rendered cube shape {cube.shape}; expected (nx, ny, nf, 2).")
    # Convert to (ny, nx, nf) for panel plotting and comparisons.
    ti = np.transpose(cube[:, :, :, 0], (1, 0, 2))
    tv = np.transpose(cube[:, :, :, 1], (1, 0, 2))
    return ti, tv, freqs, xc, yc, nx, ny


def parse_args():
    p = argparse.ArgumentParser(description="Render two models and compare output maps.")
    p.add_argument("--model1", type=Path, required=True)
    p.add_argument("--model2", type=Path, required=True)
    p.add_argument("--model1-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument("--model2-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument("--label1", type=str, default="model1")
    p.add_argument("--label2", type=str, default="model2")
    p.add_argument("--ebtel-path", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_TMP_DIR / "gximagecomputing_compare_outputs")
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument("--nx", type=int, default=None, help="Optional common output width (pixels) for both models.")
    p.add_argument("--ny", type=int, default=None, help="Optional common output height (pixels) for both models.")
    p.add_argument("--pixel-scale-arcsec", type=float, default=2.0)
    p.add_argument("--omp-threads", type=int, default=8)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if (args.nx is None) ^ (args.ny is None):
        raise ValueError("Please provide both --nx and --ny, or neither.")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    h51 = _run_render(
        model_path=args.model1,
        model_format=args.model1_format,
        output_dir=args.output_dir,
        output_name="result_model1.h5",
        ebtel_path=args.ebtel_path,
        xc=args.xc,
        yc=args.yc,
        nx=args.nx,
        ny=args.ny,
        pixel_scale_arcsec=args.pixel_scale_arcsec,
        omp_threads=args.omp_threads,
    )
    ti1, tv1, freqs, xc1, yc1, nx1, ny1 = _load_render_h5(h51)
    xc2 = args.xc if args.xc is not None else xc1
    yc2 = args.yc if args.yc is not None else yc1
    nx2 = args.nx if args.nx is not None else nx1
    ny2 = args.ny if args.ny is not None else ny1

    h52 = _run_render(
        model_path=args.model2,
        model_format=args.model2_format,
        output_dir=args.output_dir,
        output_name="result_model2.h5",
        ebtel_path=args.ebtel_path,
        xc=xc2,
        yc=yc2,
        nx=nx2,
        ny=ny2,
        pixel_scale_arcsec=args.pixel_scale_arcsec,
        omp_threads=args.omp_threads,
    )

    ti2, tv2, freqs2, _, _, _, _ = _load_render_h5(h52)
    if freqs.shape != freqs2.shape or not np.allclose(freqs, freqs2, atol=1e-6, rtol=0):
        raise ValueError(f"Frequency mismatch between model1 and model2 outputs: {freqs} vs {freqs2}")

    if ti1.shape != ti2.shape or tv1.shape != tv2.shape:
        raise ValueError(f"Shape mismatch: model1 TI/TV {ti1.shape}/{tv1.shape}, model2 TI/TV {ti2.shape}/{tv2.shape}")

    ti_abs = np.abs(ti1 - ti2)
    tv_abs = np.abs(tv1 - tv2)
    ti_pct = _percent_diff_vs_model2(ti1, ti2)
    tv_pct = _percent_diff_vs_model2(tv1, tv2)

    _panel_4x4(ti1, freqs, args.output_dir / "panel_ti_model1_4x4.png", f"TI {args.label1}")
    _panel_4x4(ti2, freqs, args.output_dir / "panel_ti_model2_4x4.png", f"TI {args.label2}")
    _panel_4x4(ti_pct, freqs, args.output_dir / "panel_ti_pct_diff_model1_vs_model2_4x4.png", "TI %DIFF")

    _panel_4x4(tv1, freqs, args.output_dir / "panel_tv_model1_4x4.png", f"TV {args.label1}")
    _panel_4x4(tv2, freqs, args.output_dir / "panel_tv_model2_4x4.png", f"TV {args.label2}")
    _panel_4x4(tv_pct, freqs, args.output_dir / "panel_tv_pct_diff_model1_vs_model2_4x4.png", "TV %DIFF")

    summary = {
        "inputs": {
            "model1": str(args.model1),
            "model2": str(args.model2),
            "label1": args.label1,
            "label2": args.label2,
        },
        "TI_model1": _stats(ti1),
        "TI_model2": _stats(ti2),
        "TI_abs_diff": _stats(ti_abs),
        "TI_pct_diff_vs_model2": _stats(ti_pct),
        "TV_model1": _stats(tv1),
        "TV_model2": _stats(tv2),
        "TV_abs_diff": _stats(tv_abs),
        "TV_pct_diff_vs_model2": _stats(tv_pct),
    }

    per_frequency = []
    for i, fghz in enumerate(freqs):
        per_frequency.append(
            {
                "freq_GHz": float(fghz),
                "TI_abs_diff": _stats(ti_abs[:, :, i]),
                "TI_pct_diff_vs_model2": _stats(ti_pct[:, :, i]),
                "TV_abs_diff": _stats(tv_abs[:, :, i]),
                "TV_pct_diff_vs_model2": _stats(tv_pct[:, :, i]),
            }
        )

    report = {"summary": summary, "per_frequency": per_frequency}
    out_json = args.output_dir / "comparison_render_outputs.json"
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("Outputs:")
    print(f"- result_model1_h5: {h51}")
    print(f"- result_model2_h5: {h52}")
    print(f"- comparison_json: {out_json}")
    print(f"- panel_ti_model1: {args.output_dir / 'panel_ti_model1_4x4.png'}")
    print(f"- panel_ti_model2: {args.output_dir / 'panel_ti_model2_4x4.png'}")
    print(f"- panel_ti_pct_diff: {args.output_dir / 'panel_ti_pct_diff_model1_vs_model2_4x4.png'}")
    print(f"- panel_tv_model1: {args.output_dir / 'panel_tv_model1_4x4.png'}")
    print(f"- panel_tv_model2: {args.output_dir / 'panel_tv_model2_4x4.png'}")
    print(f"- panel_tv_pct_diff: {args.output_dir / 'panel_tv_pct_diff_model1_vs_model2_4x4.png'}")
    print("TI abs diff stats:", summary["TI_abs_diff"])
    print("TV abs diff stats:", summary["TV_abs_diff"])


if __name__ == "__main__":
    main()
