#!/usr/bin/env python3

"""
Compare internal gximagecomputing model representations derived from:
1) pyAMPP CHR HDF5
2) GX Simulator IDL CHR SAV

This helps isolate model-conversion differences before rendering.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

import numpy as np

from gximagecomputing.radio import GXRadioImageComputing


def _is_numeric(arr: np.ndarray) -> bool:
    return np.issubdtype(arr.dtype, np.number)


def _stats_numeric(h5_arr: np.ndarray, sav_arr: np.ndarray, rel_eps: float) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    if h5_arr.shape != sav_arr.shape:
        out["shape_match"] = False
        return out

    out["shape_match"] = True
    out["exact_equal"] = bool(np.array_equal(h5_arr, sav_arr))

    finite_mask = np.isfinite(h5_arr) & np.isfinite(sav_arr)
    out["n_total"] = int(h5_arr.size)
    out["n_finite_both"] = int(np.count_nonzero(finite_mask))
    out["n_nan_h5"] = int(np.count_nonzero(np.isnan(h5_arr)))
    out["n_nan_sav"] = int(np.count_nonzero(np.isnan(sav_arr)))

    if out["n_finite_both"] == 0:
        out["allclose"] = False
        return out

    h = h5_arr[finite_mask].astype(np.float64, copy=False)
    s = sav_arr[finite_mask].astype(np.float64, copy=False)
    d = h - s
    ad = np.abs(d)
    denom = np.abs(s)
    rel_mask = denom > rel_eps
    rel = np.empty_like(ad)
    rel[:] = np.nan
    rel[rel_mask] = ad[rel_mask] / denom[rel_mask]

    out["mean_h5"] = float(np.mean(h))
    out["mean_sav"] = float(np.mean(s))
    out["mean_abs_diff"] = float(np.mean(ad))
    out["std_abs_diff"] = float(np.std(ad))
    out["max_abs_diff"] = float(np.max(ad))
    out["p95_abs_diff"] = float(np.percentile(ad, 95))

    valid_rel = rel[np.isfinite(rel)]
    if valid_rel.size > 0:
        out["mean_rel_diff"] = float(np.mean(valid_rel))
        out["max_rel_diff"] = float(np.max(valid_rel))
        out["p95_rel_diff"] = float(np.percentile(valid_rel, 95))
    else:
        out["mean_rel_diff"] = None
        out["max_rel_diff"] = None
        out["p95_rel_diff"] = None

    out["allclose"] = bool(np.allclose(h5_arr, sav_arr, rtol=1e-5, atol=1e-6, equal_nan=True))
    return out


def compare_models(
    model_h5: np.ndarray,
    model_sav: np.ndarray,
    rel_eps: float = 1e-12,
) -> Dict[str, Any]:
    fields_h5 = set(model_h5.dtype.names or [])
    fields_sav = set(model_sav.dtype.names or [])
    common = sorted(fields_h5 & fields_sav)

    report: Dict[str, Any] = {
        "summary": {},
        "fields_only_h5": sorted(fields_h5 - fields_sav),
        "fields_only_sav": sorted(fields_sav - fields_h5),
        "fields": {},
    }

    shape_mismatch = 0
    dtype_mismatch = 0
    non_numeric = 0
    not_allclose = 0

    for name in common:
        h5_val = np.asarray(model_h5[name][0])
        sav_val = np.asarray(model_sav[name][0])

        field_info: Dict[str, Any] = {
            "h5_shape": list(h5_val.shape),
            "sav_shape": list(sav_val.shape),
            "h5_dtype": str(h5_val.dtype),
            "sav_dtype": str(sav_val.dtype),
        }

        if h5_val.shape != sav_val.shape:
            shape_mismatch += 1
        if str(h5_val.dtype) != str(sav_val.dtype):
            dtype_mismatch += 1

        if _is_numeric(h5_val) and _is_numeric(sav_val):
            field_info.update(_stats_numeric(h5_val, sav_val, rel_eps=rel_eps))
            if field_info.get("allclose") is False:
                not_allclose += 1
        else:
            non_numeric += 1
            field_info["shape_match"] = bool(h5_val.shape == sav_val.shape)
            field_info["exact_equal"] = bool(np.array_equal(h5_val, sav_val))
            field_info["allclose"] = field_info["exact_equal"]
            if field_info["allclose"] is False:
                not_allclose += 1

        report["fields"][name] = field_info

    report["summary"] = {
        "n_common_fields": len(common),
        "n_fields_only_h5": len(report["fields_only_h5"]),
        "n_fields_only_sav": len(report["fields_only_sav"]),
        "n_shape_mismatch": shape_mismatch,
        "n_dtype_mismatch": dtype_mismatch,
        "n_non_numeric": non_numeric,
        "n_not_allclose": not_allclose,
    }

    # Rank biggest differences by absolute and relative metrics.
    sortable = []
    for name, info in report["fields"].items():
        sortable.append(
            (
                name,
                info.get("max_abs_diff", -1.0) if info.get("max_abs_diff") is not None else -1.0,
                info.get("max_rel_diff", -1.0) if info.get("max_rel_diff") is not None else -1.0,
            )
        )
    sortable_abs = sorted(sortable, key=lambda x: x[1], reverse=True)
    sortable_rel = sorted(sortable, key=lambda x: x[2], reverse=True)
    report["top_max_abs_diff"] = [{"field": k, "max_abs_diff": v} for k, v, _ in sortable_abs[:10]]
    report["top_max_rel_diff"] = [{"field": k, "max_rel_diff": v} for k, _, v in sortable_rel[:10]]

    return report


def parse_args():
    parser = argparse.ArgumentParser(description="Compare gximagecomputing models loaded from HDF5 vs SAV.")
    parser.add_argument(
        "--h5-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.h5",
        help="Path to CHR HDF5 file.",
    )
    parser.add_argument(
        "--sav-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.sav",
        help="Path to CHR SAV file.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Output JSON report path. Default: <system-temp>/gximage_model_compare_<timestamp>.json",
    )
    parser.add_argument(
        "--rel-eps",
        type=float,
        default=1e-12,
        help="Small epsilon for relative-difference denominator guard.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_json is None:
        args.output_json = Path(tempfile.gettempdir()) / f"gximage_model_compare_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    args.output_json.parent.mkdir(parents=True, exist_ok=True)

    # Avoid constructor to compare models without requiring native library availability.
    gxi = GXRadioImageComputing.__new__(GXRadioImageComputing)
    model_h5, _ = gxi.load_model_hdf(str(args.h5_path))
    model_sav, _ = gxi.load_model_sav(str(args.sav_path))

    report = compare_models(model_h5, model_sav, rel_eps=args.rel_eps)
    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("Outputs:")
    print(f"- comparison_json: {args.output_json}")
    print("Summary:")
    for k, v in report["summary"].items():
        print(f"- {k}: {v}")
    print("Top max_abs_diff fields:")
    for row in report["top_max_abs_diff"][:5]:
        print(f"- {row['field']}: {row['max_abs_diff']}")
    print("Top max_rel_diff fields:")
    for row in report["top_max_rel_diff"][:5]:
        print(f"- {row['field']}: {row['max_rel_diff']}")


if __name__ == "__main__":
    main()
