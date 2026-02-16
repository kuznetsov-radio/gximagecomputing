#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict

import numpy as np

from gximagecomputing.radio import GXRadioImageComputing


FIELDS = [
    "dz",
    "Bx",
    "By",
    "Bz",
    "chromo_n0",
    "chromo_np",
    "chromo_nHI",
    "chromo_T0",
    "corona_Bavg",
    "corona_L",
    "chromo_uniform_Bavg",
    "chromo_uniform_L",
    "corona_ID1",
    "corona_ID2",
    "chromo_uniform_ID1",
    "chromo_uniform_ID2",
]


def _stats(h5: np.ndarray, sav: np.ndarray, eps: float = 1e-12) -> Dict[str, float]:
    diff = h5 - sav
    abs_diff = np.abs(diff)
    denom = np.abs(sav)
    rel = np.full(abs_diff.shape, np.nan, dtype=np.float64)
    m = denom > eps
    rel[m] = abs_diff[m] / denom[m]
    return {
        "mean_h5": float(np.nanmean(h5)),
        "mean_sav": float(np.nanmean(sav)),
        "mean_abs_diff": float(np.nanmean(abs_diff)),
        "std_abs_diff": float(np.nanstd(abs_diff)),
        "max_abs_diff": float(np.nanmax(abs_diff)),
        "mean_rel_diff": float(np.nanmean(rel)),
        "p95_rel_diff": float(np.nanpercentile(rel, 95)),
        "max_rel_diff": float(np.nanmax(rel)),
    }


def parse_args():
    p = argparse.ArgumentParser(description="Compare normalized render-input arrays (H5 vs SAV).")
    p.add_argument(
        "--h5-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.h5",
    )
    p.add_argument(
        "--sav-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.sav",
    )
    p.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Default: <system-temp>/gximage_render_inputs_compare_<timestamp>.json",
    )
    return p.parse_args()


def main():
    args = parse_args()
    if args.output_json is None:
        args.output_json = Path(tempfile.gettempdir()) / f"gximage_render_inputs_compare_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

    g = GXRadioImageComputing.__new__(GXRadioImageComputing)
    mh5, _ = g.load_model_hdf(str(args.h5_path))
    msav, _ = g.load_model_sav(str(args.sav_path))

    report = {"fields": {}, "summary": {}}
    n_not_close = 0

    for f in FIELDS:
        h = np.asarray(mh5[f][0])
        s = np.asarray(msav[f][0])
        entry = {
            "shape_h5": list(h.shape),
            "shape_sav": list(s.shape),
            "dtype_h5": str(h.dtype),
            "dtype_sav": str(s.dtype),
        }
        if h.shape != s.shape:
            entry["shape_match"] = False
            report["fields"][f] = entry
            n_not_close += 1
            continue
        entry["shape_match"] = True
        if np.issubdtype(h.dtype, np.number) and np.issubdtype(s.dtype, np.number):
            entry.update(_stats(h.astype(np.float64), s.astype(np.float64)))
            entry["allclose"] = bool(np.allclose(h, s, rtol=1e-5, atol=1e-6, equal_nan=True))
        else:
            entry["allclose"] = bool(np.array_equal(h, s))
        if not entry["allclose"]:
            n_not_close += 1
        report["fields"][f] = entry

    ranked_abs = sorted(
        (
            (k, v.get("max_abs_diff", -1.0))
            for k, v in report["fields"].items()
        ),
        key=lambda x: x[1],
        reverse=True,
    )
    ranked_rel = sorted(
        (
            (k, v.get("max_rel_diff", -1.0))
            for k, v in report["fields"].items()
        ),
        key=lambda x: x[1],
        reverse=True,
    )

    report["summary"] = {
        "n_fields": len(FIELDS),
        "n_not_allclose": n_not_close,
        "top_max_abs_diff": ranked_abs[:6],
        "top_max_rel_diff": ranked_rel[:6],
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, "w", encoding="utf-8") as fp:
        json.dump(report, fp, indent=2)

    print("Outputs:")
    print(f"- comparison_json: {args.output_json}")
    print("Summary:")
    for k, v in report["summary"].items():
        print(f"- {k}: {v}")


if __name__ == "__main__":
    main()
