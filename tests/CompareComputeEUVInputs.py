#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import readsav


TOP_VARS = ("model", "ebtel", "response", "simbox", "coronaparms", "outspace", "shtable")


def _parse_args() -> argparse.Namespace:
    base = Path(tempfile.gettempdir()) / "gximagecomputing_validation_groundtruth"
    p = argparse.ArgumentParser(description="Compare Python-vs-IDL ComputeEUV inputs (pre-DLL-call dumps).")
    p.add_argument("--python-dir", type=Path, default=base / "computeeuv_inputs_python")
    p.add_argument("--idl-sav", type=Path, default=base / "computeeuv_inputs_idl.sav")
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-9)
    p.add_argument(
        "--output-json",
        type=Path,
        default=base / f"computeeuv_inputs_compare_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
    )
    return p.parse_args()


def _decode_text_scalar(v: Any) -> str:
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace")
    return str(v)


def _unwrap_singleton_object(v: Any) -> Any:
    x = v
    for _ in range(8):
        if isinstance(x, np.ndarray) and x.dtype == object and x.size == 1:
            x = x.reshape(-1)[0]
            continue
        break
    return x


def _is_structured(v: Any) -> bool:
    if isinstance(v, np.ndarray):
        return v.dtype.names is not None
    if isinstance(v, (np.void, np.record)):
        return getattr(v, "dtype", None) is not None and v.dtype.names is not None
    return False


def _as_structured_array(v: Any) -> np.ndarray:
    if isinstance(v, np.ndarray) and v.dtype.names is not None:
        return v
    if isinstance(v, (np.void, np.record)) and v.dtype.names is not None:
        return np.asarray(v).reshape(())
    raise TypeError(f"Not a structured value: {type(v)}")


def _field_map(dtype_names: tuple[str, ...] | None) -> dict[str, str]:
    return {str(n).lower(): str(n) for n in (dtype_names or ())}


def _numeric_stats(a: np.ndarray) -> dict[str, Any]:
    x = np.asarray(a)
    out: dict[str, Any] = {
        "shape": list(x.shape),
        "dtype": str(x.dtype),
        "size": int(x.size),
    }
    if x.size == 0:
        return out
    if np.issubdtype(x.dtype, np.number):
        xf = x.astype(np.float64, copy=False)
        finite = np.isfinite(xf)
        out.update(
            {
                "finite": int(finite.sum()),
                "nan": int(np.isnan(xf).sum()),
                "inf": int(np.isinf(xf).sum()),
                "min": float(np.nanmin(xf)),
                "max": float(np.nanmax(xf)),
                "mean": float(np.nanmean(xf)),
            }
        )
    return out


def _compare_leaf(py_v: Any, idl_v: Any, rtol: float, atol: float) -> dict[str, Any]:
    py_v = _unwrap_singleton_object(py_v)
    idl_v = _unwrap_singleton_object(idl_v)

    py_a = np.asarray(py_v)
    idl_a = np.asarray(idl_v)
    report: dict[str, Any] = {
        "kind": "leaf",
        "python": _numeric_stats(py_a),
        "idl": _numeric_stats(idl_a),
        "shape_match": list(py_a.shape) == list(idl_a.shape),
        "dtype_python": str(py_a.dtype),
        "dtype_idl": str(idl_a.dtype),
    }

    if py_a.dtype.names is not None or idl_a.dtype.names is not None:
        report["kind"] = "unexpected_struct_leaf"
        report["allclose"] = False
        return report

    if py_a.shape != idl_a.shape:
        report["allclose"] = False
        if py_a.ndim == idl_a.ndim and tuple(py_a.shape) == tuple(reversed(idl_a.shape)):
            report["hint"] = "shapes are exact reverse-order match"
        return report

    if np.issubdtype(py_a.dtype, np.number) and np.issubdtype(idl_a.dtype, np.number):
        py_n = py_a.astype(np.float64, copy=False)
        idl_n = idl_a.astype(np.float64, copy=False)
        same = np.allclose(py_n, idl_n, rtol=rtol, atol=atol, equal_nan=True)
        diff = py_n - idl_n
        abs_diff = np.abs(diff)
        finite_both = np.isfinite(py_n) & np.isfinite(idl_n)
        report["allclose"] = bool(same)
        report["numeric"] = {
            "rtol": rtol,
            "atol": atol,
            "mean_abs_diff": float(np.nanmean(abs_diff)),
            "max_abs_diff": float(np.nanmax(abs_diff)),
            "mean_signed_diff": float(np.nanmean(diff)),
            "finite_overlap": int(finite_both.sum()),
        }
        return report

    if py_a.dtype.kind in {"S", "U", "O"} or idl_a.dtype.kind in {"S", "U", "O"}:
        py_txt = [_decode_text_scalar(x) for x in py_a.reshape(-1)]
        idl_txt = [_decode_text_scalar(x) for x in idl_a.reshape(-1)]
        report["allclose"] = bool(py_txt == idl_txt)
        if not report["allclose"]:
            report["python_values"] = py_txt[:20]
            report["idl_values"] = idl_txt[:20]
        return report

    report["allclose"] = bool(np.array_equal(py_a, idl_a))
    return report


def _compare_value(py_v: Any, idl_v: Any, rtol: float, atol: float) -> dict[str, Any]:
    py_v = _unwrap_singleton_object(py_v)
    idl_v = _unwrap_singleton_object(idl_v)

    py_is_struct = _is_structured(py_v)
    idl_is_struct = _is_structured(idl_v)
    if py_is_struct != idl_is_struct:
        return {
            "kind": "type_mismatch",
            "python_type": str(type(py_v)),
            "idl_type": str(type(idl_v)),
            "allclose": False,
        }
    if not py_is_struct:
        return _compare_leaf(py_v, idl_v, rtol=rtol, atol=atol)

    py_s = _as_structured_array(py_v)
    idl_s = _as_structured_array(idl_v)
    report: dict[str, Any] = {
        "kind": "struct",
        "python_shape": list(py_s.shape),
        "idl_shape": list(idl_s.shape),
        "shape_match": list(py_s.shape) == list(idl_s.shape),
        "python_dtype_names": list(py_s.dtype.names or ()),
        "idl_dtype_names": list(idl_s.dtype.names or ()),
        "fields": {},
        "allclose": True,
    }
    if py_s.shape != idl_s.shape:
        report["allclose"] = False

    py_fields = _field_map(py_s.dtype.names)
    idl_fields = _field_map(idl_s.dtype.names)
    missing_in_idl = sorted(k for k in py_fields if k not in idl_fields)
    missing_in_py = sorted(k for k in idl_fields if k not in py_fields)
    if missing_in_idl:
        report["missing_fields_in_idl"] = [py_fields[k] for k in missing_in_idl]
        report["allclose"] = False
    if missing_in_py:
        report["missing_fields_in_python"] = [idl_fields[k] for k in missing_in_py]
        report["allclose"] = False

    for key in sorted(set(py_fields) & set(idl_fields)):
        py_name = py_fields[key]
        idl_name = idl_fields[key]
        sub = _compare_value(py_s[py_name], idl_s[idl_name], rtol=rtol, atol=atol)
        report["fields"][py_name] = {"idl_field": idl_name, **sub}
        if not bool(sub.get("allclose", False)):
            report["allclose"] = False
    return report


def _load_python_dir(py_dir: Path) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in TOP_VARS:
        fp = py_dir / f"{name}.npy"
        if fp.exists():
            out[name] = np.load(fp, allow_pickle=False)
    manifest = py_dir / "manifest.json"
    if manifest.exists():
        out["_manifest"] = json.loads(manifest.read_text(encoding="utf-8"))
    return out


def _load_idl_sav(idl_sav: Path) -> dict[str, Any]:
    return readsav(str(idl_sav), python_dict=True, verbose=False)


def _count_failures(node: Any) -> int:
    if not isinstance(node, dict):
        return 0
    n = 0
    if "allclose" in node and node.get("allclose") is False:
        n += 1
    for v in node.values():
        if isinstance(v, dict):
            n += _count_failures(v)
    return n


def main() -> None:
    args = _parse_args()
    py = _load_python_dir(args.python_dir)
    idl = _load_idl_sav(args.idl_sav)

    report: dict[str, Any] = {
        "inputs": {"python_dir": str(args.python_dir), "idl_sav": str(args.idl_sav)},
        "top": {},
        "summary": {},
    }
    if "_manifest" in py:
        report["python_manifest"] = py["_manifest"]

    top_fail = 0
    missing_py = []
    missing_idl = []
    for name in TOP_VARS:
        py_v = py.get(name, None)
        idl_key = next((k for k in idl.keys() if str(k).lower() == name.lower()), None)
        if py_v is None:
            missing_py.append(name)
            continue
        if idl_key is None:
            missing_idl.append(name)
            continue
        cmp_node = _compare_value(py_v, idl[idl_key], rtol=args.rtol, atol=args.atol)
        report["top"][name] = cmp_node
        if not bool(cmp_node.get("allclose", False)):
            top_fail += 1

    report["summary"] = {
        "missing_top_vars_in_python_dump": missing_py,
        "missing_top_vars_in_idl_sav": missing_idl,
        "n_top_compared": len(report["top"]),
        "n_top_not_allclose": top_fail,
        "n_recursive_fail_nodes": _count_failures(report),
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, "w", encoding="utf-8") as fp:
        json.dump(report, fp, indent=2)

    print("ComputeEUV input comparison:")
    print(f"- python_dir: {args.python_dir}")
    print(f"- idl_sav: {args.idl_sav}")
    print(f"- output_json: {args.output_json}")
    print("Top-level results:")
    for name in TOP_VARS:
        if name in report["top"]:
            print(f"- {name}: allclose={report['top'][name].get('allclose')}")
    if missing_py:
        print(f"- missing_in_python: {missing_py}")
    if missing_idl:
        print(f"- missing_in_idl: {missing_idl}")
    print(f"- n_top_not_allclose: {report['summary']['n_top_not_allclose']}")


if __name__ == "__main__":
    main()
