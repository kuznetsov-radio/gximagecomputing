#!/usr/bin/env python3

"""
Regression guard for SAV <-> H5 model parity.

Usage examples:
  python tests/RegressionRoundTripSavH5.py
  python tests/RegressionRoundTripSavH5.py --rebuild-h5
  python tests/RegressionRoundTripSavH5.py --sav-path /path/to/model.sav --h5-path /path/to/model.h5
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

# Avoid SunPy config write failures in restricted environments.
os.environ.setdefault(
    "SUNPY_CONFIGDIR",
    str(Path(tempfile.gettempdir()) / "gximagecomputing_sunpy_config"),
)

from gximagecomputing.io.sav_to_h5 import build_h5_from_sav
from gximagecomputing.radio import GXRadioImageComputing


def _parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    test_data = root / "test_data"
    p = argparse.ArgumentParser(description="Strict regression check for SAV/H5 parity.")
    p.add_argument(
        "--sav-path",
        type=Path,
        default=test_data / "test.chr.sav",
        help="Reference SAV file.",
    )
    p.add_argument(
        "--h5-path",
        type=Path,
        default=test_data / "test.chr.h5",
        help="H5 file to compare against SAV.",
    )
    p.add_argument(
        "--rebuild-h5",
        action="store_true",
        help="Rebuild --h5-path from --sav-path before comparison.",
    )
    p.add_argument(
        "--atol",
        type=float,
        default=0.0,
        help="Absolute tolerance for float fields.",
    )
    p.add_argument(
        "--rtol",
        type=float,
        default=0.0,
        help="Relative tolerance for float fields.",
    )
    return p.parse_args()


def _compare_models(model_h5: np.ndarray, model_sav: np.ndarray, atol: float, rtol: float) -> Dict[str, str]:
    failures: Dict[str, str] = {}

    names_h5 = set(model_h5.dtype.names or [])
    names_sav = set(model_sav.dtype.names or [])
    if names_h5 != names_sav:
        failures["dtype_names"] = (
            f"Field sets differ. only_h5={sorted(names_h5 - names_sav)}, "
            f"only_sav={sorted(names_sav - names_h5)}"
        )
        return failures

    for name in sorted(names_h5):
        h = np.asarray(model_h5[name][0])
        s = np.asarray(model_sav[name][0])

        if h.shape != s.shape:
            failures[name] = f"shape mismatch: h5={h.shape}, sav={s.shape}"
            continue

        if np.issubdtype(h.dtype, np.number) and np.issubdtype(s.dtype, np.number):
            if np.issubdtype(h.dtype, np.floating) or np.issubdtype(s.dtype, np.floating):
                if not np.allclose(h, s, atol=atol, rtol=rtol, equal_nan=True):
                    max_abs = float(np.max(np.abs(h.astype(np.float64) - s.astype(np.float64))))
                    failures[name] = f"float mismatch: max_abs={max_abs:.6g}, atol={atol}, rtol={rtol}"
            else:
                if not np.array_equal(h, s):
                    failures[name] = "integer mismatch"
        else:
            if not np.array_equal(h, s):
                failures[name] = "non-numeric mismatch"

    return failures


def main() -> int:
    args = _parse_args()
    sav_path = args.sav_path.expanduser().resolve()
    h5_path = args.h5_path.expanduser().resolve()

    if not sav_path.exists():
        print(f"ERROR: missing SAV file: {sav_path}", file=sys.stderr)
        return 2

    if args.rebuild_h5:
        h5_path.parent.mkdir(parents=True, exist_ok=True)
        build_h5_from_sav(sav_path=sav_path, out_h5=h5_path, template_h5=None)
        print(f"Rebuilt H5: {h5_path}")

    if not h5_path.exists():
        print(f"ERROR: missing H5 file: {h5_path}", file=sys.stderr)
        return 2

    # CHR-only regression: reject non-CHR outputs even though converter is stage-agnostic.
    import h5py
    with h5py.File(h5_path, "r") as f:
        if "chromo" not in f:
            print(
                f"ERROR: {h5_path} is not CHR-compatible (missing /chromo). "
                "Use this regression only with CHR models.",
                file=sys.stderr,
            )
            return 2
        stage = None
        if "metadata" in f and "stage" in f["metadata"]:
            raw = f["metadata"]["stage"][()]
            if isinstance(raw, bytes):
                stage = raw.decode("utf-8", "ignore").upper()
            else:
                stage = str(raw).upper()
        if stage is not None and stage != "CHR":
            print(
                f"ERROR: {h5_path} has metadata stage={stage}, expected CHR for this regression.",
                file=sys.stderr,
            )
            return 2

    # Avoid constructor so the native renderer library is not required for this regression check.
    gxi = GXRadioImageComputing.__new__(GXRadioImageComputing)
    model_h5, _ = gxi.load_model_hdf(str(h5_path))
    model_sav, _ = gxi.load_model_sav(str(sav_path))

    failures = _compare_models(model_h5, model_sav, atol=args.atol, rtol=args.rtol)
    if failures:
        print("FAILED parity regression:")
        for k, v in failures.items():
            print(f"- {k}: {v}")
        return 1

    print("PASS parity regression: SAV and H5 models are equivalent.")
    print(f"- sav: {sav_path}")
    print(f"- h5: {h5_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
