#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


def _default_idl_file() -> Path:
    # Canonical source file path from the user request (fixed typo: M_720s).
    return Path(
        "/Users/gelu/Library/CloudStorage/Dropbox/@Projects/sim4fasr/gx_models/2025-11-26/"
        "hmi.M_720s.20251126_153431.W28S12CR.CEA.NAS.CHR.sav"
    )


def _default_test_data_dir() -> Path:
    return Path(
        "/Users/gelu/Library/CloudStorage/Dropbox/@Projects/@SUNCAST-ORG/gximagecomputing/test_data"
    )


def _import_builder(repo_root: Path):
    src_dir = repo_root / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))
    from gximagecomputing.io.sav_to_h5 import build_h5_from_sav

    return build_h5_from_sav


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Seed gximagecomputing/test_data from an IDL CHR SAV file by copying to "
            "test.chr.sav and converting to test.chr.h5."
        )
    )
    p.add_argument(
        "--idl-file",
        type=Path,
        default=_default_idl_file(),
        help="Input CHR SAV file path.",
    )
    p.add_argument(
        "--test-data-dir",
        type=Path,
        default=_default_test_data_dir(),
        help="Destination test_data directory.",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    build_h5_from_sav = _import_builder(repo_root)

    idl_file = args.idl_file.expanduser().resolve()
    test_data_dir = args.test_data_dir.expanduser().resolve()
    test_data_dir.mkdir(parents=True, exist_ok=True)

    if not idl_file.exists():
        raise FileNotFoundError(f"Input SAV not found: {idl_file}")

    out_sav = test_data_dir / "test.chr.sav"
    out_h5 = test_data_dir / "test.chr.h5"

    shutil.copy2(idl_file, out_sav)
    build_h5_from_sav(sav_path=out_sav, out_h5=out_h5, template_h5=None)

    print("Created test data:")
    print(f"- {out_sav}")
    print(f"- {out_h5}")


if __name__ == "__main__":
    main()
