#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy.io import readsav


def _decode(value: Any) -> Any:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", "ignore")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray) and value.shape == ():
        return _decode(value.item())
    return value


def _load_idl_dump(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    payload = readsav(str(path), python_dict=True, verbose=False)
    meta = {
        "source_file": _decode(payload["source_file"]),
        "source_name": _decode(payload["source_name"]),
        "source_kind": _decode(payload["source_kind"]),
        "case_name": _decode(payload["case_name"]),
        "override_applied": bool(np.asarray(payload["override_applied"]).reshape(-1)[0]),
        "override_dsun_cm": float(np.asarray(payload["override_dsun_cm"]).reshape(-1)[0]),
        "override_lonc_deg": float(np.asarray(payload["override_lonc_deg"]).reshape(-1)[0]),
        "override_b0sun_deg": float(np.asarray(payload["override_b0sun_deg"]).reshape(-1)[0]),
    }
    model = np.asarray(payload["model"])
    name_map = {
        "NX": "Nx",
        "NY": "Ny",
        "NZ": "Nz",
        "CHROMO_LAYERS": "chromo_layers",
        "CORONA_LAYERS": "corona_layers",
        "CORONA_BASE": "corona_base",
        "DSUN": "DSun",
        "RSUN": "RSun",
        "B0SUN": "b0Sun",
        "LONC": "lonC",
        "LATC": "latC",
        "DX": "dx",
        "DY": "dy",
        "DZ_UNIFORM": "dz_uniform",
        "OBSTIME": "obstime",
        "DZ": "dz",
        "BX": "Bx",
        "BY": "By",
        "BZ": "Bz",
        "CHROMO_N0": "chromo_n0",
        "CHROMO_NP": "chromo_np",
        "CHROMO_NHI": "chromo_nHI",
        "CHROMO_T0": "chromo_T0",
        "CORONA_BAVG": "corona_Bavg",
        "CORONA_L": "corona_L",
        "CHROMO_UNIFORM_BAVG": "chromo_uniform_Bavg",
        "CHROMO_UNIFORM_L": "chromo_uniform_L",
        "VOXELID": "VoxelID",
        "CORONA_ID1": "corona_ID1",
        "CORONA_ID2": "corona_ID2",
        "CHROMO_UNIFORM_ID1": "chromo_uniform_ID1",
        "CHROMO_UNIFORM_ID2": "chromo_uniform_ID2",
    }
    fields = {
        name_map.get(name, name): np.asarray(model[name][0]) for name in (model.dtype.names or ())
    }
    return meta, fields


def _load_python_dump(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    with h5py.File(path, "r") as f:
        meta = {key: _decode(f.attrs[key]) for key in f.attrs.keys()}
        fields = {name: np.asarray(f["model"][name][...]) for name in f["model"].keys()}
    return meta, fields


def _numeric_stats(py_value: np.ndarray, idl_value: np.ndarray) -> dict[str, Any]:
    info: dict[str, Any] = {
        "shape_match": bool(py_value.shape == idl_value.shape),
        "exact_equal": bool(np.array_equal(py_value, idl_value)),
    }
    if py_value.shape != idl_value.shape:
        return info
    diff = np.asarray(py_value, dtype=np.float64) - np.asarray(idl_value, dtype=np.float64)
    abs_diff = np.abs(diff)
    info["allclose"] = bool(np.allclose(py_value, idl_value, rtol=1e-6, atol=1e-6, equal_nan=True))
    info["max_abs_diff"] = float(np.nanmax(abs_diff))
    info["mean_abs_diff"] = float(np.nanmean(abs_diff))
    return info


def compare_field_maps(python_fields: dict[str, np.ndarray], idl_fields: dict[str, np.ndarray]) -> dict[str, Any]:
    python_names = set(python_fields)
    idl_names = set(idl_fields)
    common = sorted(python_names & idl_names)
    report: dict[str, Any] = {
        "fields_only_python": sorted(python_names - idl_names),
        "fields_only_idl": sorted(idl_names - python_names),
        "fields": {},
    }

    n_shape_mismatch = 0
    n_not_allclose = 0
    for name in common:
        py_value = python_fields[name]
        idl_value = idl_fields[name]
        info: dict[str, Any] = {
            "python_shape": list(py_value.shape),
            "idl_shape": list(idl_value.shape),
            "python_dtype": str(py_value.dtype),
            "idl_dtype": str(idl_value.dtype),
        }
        if np.issubdtype(py_value.dtype, np.number) and np.issubdtype(idl_value.dtype, np.number):
            info.update(_numeric_stats(py_value, idl_value))
        else:
            info["shape_match"] = bool(py_value.shape == idl_value.shape)
            info["exact_equal"] = bool(np.array_equal(py_value, idl_value))
            info["allclose"] = info["exact_equal"]
        if not info["shape_match"]:
            n_shape_mismatch += 1
        if not info["allclose"]:
            n_not_allclose += 1
        report["fields"][name] = info

    report["summary"] = {
        "n_common_fields": len(common),
        "n_fields_only_python": len(report["fields_only_python"]),
        "n_fields_only_idl": len(report["fields_only_idl"]),
        "n_shape_mismatch": n_shape_mismatch,
        "n_not_allclose": n_not_allclose,
    }
    return report


def parse_args():
    parser = argparse.ArgumentParser(description="Compare IDL and Python loader dump outputs.")
    parser.add_argument(
        "--idl-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "gximagecomputing_loader_dumps" / "idl",
        help="Directory containing *.loaderdump.sav files from DumpIDLLoadGXmodelParity.pro.",
    )
    parser.add_argument(
        "--python-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "gximagecomputing_loader_dumps" / "python",
        help="Directory containing *.loaderdump.h5 files from DumpPythonLoadModelParity.py.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path(tempfile.gettempdir()) / f"gximagecomputing_loader_dump_parity_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        help="Output JSON report path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    idl_files = {path.name.removesuffix(".loaderdump.sav"): path for path in sorted(args.idl_dir.glob("*.loaderdump.sav"))}
    python_files = {
        path.name.removesuffix(".loaderdump.h5"): path for path in sorted(args.python_dir.glob("*.loaderdump.h5"))
    }

    pairs = sorted(set(idl_files) & set(python_files))
    report: dict[str, Any] = {
        "idl_only": sorted(set(idl_files) - set(python_files)),
        "python_only": sorted(set(python_files) - set(idl_files)),
        "pairs": {},
    }

    for key in pairs:
        idl_meta, idl_fields = _load_idl_dump(idl_files[key])
        py_meta, py_fields = _load_python_dump(python_files[key])
        pair_report = compare_field_maps(py_fields, idl_fields)
        pair_report["idl_meta"] = idl_meta
        pair_report["python_meta"] = py_meta

        if py_meta.get("override_applied"):
            expected_dsun = float(py_meta["override_dsun_cm"])
            expected_lonc = float(py_meta["override_lonc_deg"])
            expected_b0 = float(py_meta["override_b0sun_deg"])
            pair_report["override_checks"] = {
                "python_DSun_matches": bool(np.isclose(float(py_fields["DSun"]), expected_dsun)),
                "python_lonC_matches": bool(np.isclose(float(py_fields["lonC"]), expected_lonc)),
                "python_b0Sun_matches": bool(np.isclose(float(py_fields["b0Sun"]), expected_b0)),
                "idl_DSun_matches": bool(np.isclose(float(idl_fields["DSun"]), expected_dsun)),
                "idl_lonC_matches": bool(np.isclose(float(idl_fields["lonC"]), expected_lonc)),
                "idl_b0Sun_matches": bool(np.isclose(float(idl_fields["b0Sun"]), expected_b0)),
            }

        report["pairs"][key] = pair_report

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(_decode(report), indent=2), encoding="utf-8")

    print("Outputs:")
    print(f"- comparison_json: {args.output_json}")
    print("Summary:")
    print(f"- n_pairs: {len(pairs)}")
    print(f"- n_idl_only: {len(report['idl_only'])}")
    print(f"- n_python_only: {len(report['python_only'])}")
    for key in pairs:
        summary = report["pairs"][key]["summary"]
        print(f"- {key}: n_not_allclose={summary['n_not_allclose']}, n_shape_mismatch={summary['n_shape_mismatch']}")


if __name__ == "__main__":
    main()
