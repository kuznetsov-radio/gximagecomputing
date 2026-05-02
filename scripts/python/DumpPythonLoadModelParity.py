#!/usr/bin/env python3

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from gxrender.io.model import load_model_hdf_with_metadata, load_model_sav_with_metadata
from gxrender.utils.test_data import find_default_model_file


OVERRIDE_DSUN_CM = 1.4321098765e13
OVERRIDE_LONC_DEG = 23.456789
OVERRIDE_B0SUN_DEG = -6.54321


def _coerce_attr_value(value: Any):
    if hasattr(value, "isot"):
        return str(value.isot)
    if isinstance(value, (bytes, np.bytes_)):
        return np.bytes_(value)
    if isinstance(value, str):
        return np.bytes_(value)
    arr = np.asarray(value)
    if arr.shape == ():
        return arr.item()
    return arr


def _write_dump(
    out_path: Path,
    *,
    source_path: Path,
    case_name: str,
    override_applied: bool,
    recompute_observer_ephemeris: bool,
    observer_name: str | None,
    override_dsun_cm: float | None,
    override_lonc_deg: float | None,
    override_b0sun_deg: float | None,
    metadata: dict[str, Any],
    model: np.ndarray,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as f:
        f.attrs["source_file"] = np.bytes_(str(source_path))
        f.attrs["source_name"] = np.bytes_(source_path.name)
        f.attrs["source_kind"] = np.bytes_(source_path.suffix.lower().lstrip("."))
        f.attrs["case_name"] = np.bytes_(case_name)
        f.attrs["override_applied"] = int(override_applied)
        f.attrs["recompute_observer_ephemeris"] = int(recompute_observer_ephemeris)
        if observer_name is not None:
            f.attrs["observer_name"] = np.bytes_(str(observer_name))
        if override_dsun_cm is not None:
            f.attrs["override_dsun_cm"] = float(override_dsun_cm)
        if override_lonc_deg is not None:
            f.attrs["override_lonc_deg"] = float(override_lonc_deg)
        if override_b0sun_deg is not None:
            f.attrs["override_b0sun_deg"] = float(override_b0sun_deg)

        g_header = f.create_group("header")
        for key, value in sorted(metadata.items()):
            try:
                g_header.attrs[key] = _coerce_attr_value(value)
            except TypeError:
                g_header.attrs[key] = np.bytes_(str(value))

        g_model = f.create_group("model")
        for field in model.dtype.names or ():
            g_model.create_dataset(field, data=np.asarray(model[field][0]))


def _dump_case(source_path: Path, out_dir: Path, case_name: str, overrides: dict[str, float]) -> Path:
    if source_path.suffix.lower() == ".sav":
        model, _model_dt, metadata = load_model_sav_with_metadata(str(source_path), **overrides)
    elif source_path.suffix.lower() in {".h5", ".hdf5"}:
        model, _model_dt, metadata = load_model_hdf_with_metadata(str(source_path), **overrides)
    else:
        raise ValueError(f"Unsupported source file: {source_path}")

    out_path = out_dir / f"{source_path.name}.{case_name}.loaderdump.h5"
    _write_dump(
        out_path,
        source_path=source_path,
        case_name=case_name,
        override_applied=bool(
            overrides.get("DSun") is not None
            or overrides.get("lonC") is not None
            or overrides.get("b0Sun") is not None
        ),
        recompute_observer_ephemeris=bool(overrides.get("recompute_observer_ephemeris", False)),
        observer_name=overrides.get("observer_name"),
        override_dsun_cm=overrides.get("DSun"),
        override_lonc_deg=overrides.get("lonC"),
        override_b0sun_deg=overrides.get("b0Sun"),
        metadata=metadata,
        model=model,
    )
    return out_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Dump DLL-ready Python loader outputs for all *chr* models in test_data."
    )
    parser.add_argument(
        "--test-data-dir",
        type=Path,
        default=None,
        help="Directory containing *chr*.sav and *chr*.h5 inputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "gximagecomputing_loader_dumps" / "python",
        help="Directory for per-case HDF5 dump outputs.",
    )
    parser.add_argument("--override-dsun-cm", type=float, default=OVERRIDE_DSUN_CM)
    parser.add_argument("--override-lonc-deg", type=float, default=OVERRIDE_LONC_DEG)
    parser.add_argument("--override-b0sun-deg", type=float, default=OVERRIDE_B0SUN_DEG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.test_data_dir is None:
        args.test_data_dir = find_default_model_file().parent
    files = sorted(args.test_data_dir.glob("*chr*.sav")) + sorted(args.test_data_dir.glob("*chr*.h5"))
    cases = [
        ("default", {}),
        (
            "recompute_earth",
            {
                "recompute_observer_ephemeris": True,
                "observer_name": "earth",
            },
        ),
        (
            "override",
            {
                "DSun": float(args.override_dsun_cm),
                "lonC": float(args.override_lonc_deg),
                "b0Sun": float(args.override_b0sun_deg),
            },
        ),
    ]

    written: list[Path] = []
    for source_path in files:
        for case_name, overrides in cases:
            written.append(_dump_case(source_path, args.output_dir, case_name, overrides))

    print("Outputs:")
    for path in written:
        print(f"- {path}")


if __name__ == "__main__":
    main()
