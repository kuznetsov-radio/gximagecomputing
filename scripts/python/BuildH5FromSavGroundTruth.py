#!/usr/bin/env python3

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

from gxrender.io.sav_to_h5 import build_h5_from_sav
from gxrender.utils.test_data import find_default_model_file


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Legacy test wrapper around gxrender.io.sav_to_h5.build_h5_from_sav."
    )
    p.add_argument(
        "--template-h5",
        type=Path,
        default=None,
        help="Optional template H5 copied before writing output.",
    )
    p.add_argument(
        "--sav-path",
        type=Path,
        default=None,
    )
    p.add_argument(
        "--out-h5",
        type=Path,
        default=Path(tempfile.gettempdir()) / "chr_from_sav_groundtruth.h5",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.sav_path is None:
        args.sav_path = find_default_model_file(".sav")
    template_h5 = args.template_h5.expanduser().resolve() if args.template_h5 else None
    out_h5 = build_h5_from_sav(
        sav_path=args.sav_path,
        out_h5=args.out_h5,
        template_h5=template_h5 if template_h5 and template_h5.exists() else None,
    )
    print("Outputs:")
    print(f"- out_h5: {out_h5}")


if __name__ == "__main__":
    main()
