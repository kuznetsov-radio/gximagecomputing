#!/usr/bin/env python3

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from gximagecomputing.io.sav_to_h5 import build_h5_from_sav


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Legacy test wrapper around gximagecomputing.io.sav_to_h5.build_h5_from_sav."
    )
    p.add_argument(
        "--template-h5",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.h5",
        help="Optional template H5 copied before writing output.",
    )
    p.add_argument(
        "--sav-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.sav",
    )
    p.add_argument(
        "--out-h5",
        type=Path,
        default=Path(tempfile.gettempdir()) / "chr_from_sav_groundtruth.h5",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_h5 = build_h5_from_sav(
        sav_path=args.sav_path,
        out_h5=args.out_h5,
        template_h5=args.template_h5 if args.template_h5 and args.template_h5.exists() else None,
    )
    print("Outputs:")
    print(f"- out_h5: {out_h5}")


if __name__ == "__main__":
    main()
