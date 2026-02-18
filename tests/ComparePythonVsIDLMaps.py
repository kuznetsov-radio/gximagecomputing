#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import tempfile
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav


DEFAULT_OUTDIR = Path(tempfile.gettempdir()) / "gximagecomputing_validation_groundtruth"
DEFAULT_PY_H5 = DEFAULT_OUTDIR / "test.chr.sav_py_mw_maps.h5"
DEFAULT_IDL_SAV = DEFAULT_OUTDIR / "idl_mw_maps.sav"


def _percent_diff(test: np.ndarray, truth: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    out = np.zeros_like(test, dtype=np.float64)
    denom = np.abs(truth)
    mask = denom > eps
    out[mask] = 100.0 * (test[mask] - truth[mask]) / denom[mask]
    return out


def _sym_diff_pm1(test: np.ndarray, truth: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    # Bounded symmetric difference in [-1, 1] without masking.
    return (test - truth) / (np.abs(test) + np.abs(truth) + eps)


def _sum_diff(test: np.ndarray, truth: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    # User-requested form; not bounded when test+truth is near zero.
    return (test - truth) / (test + truth + eps)


def _stats(a: np.ndarray) -> dict[str, float]:
    x = np.asarray(a, dtype=np.float64)
    return {
        "mean": float(np.nanmean(x)),
        "std": float(np.nanstd(x)),
        "min": float(np.nanmin(x)),
        "max": float(np.nanmax(x)),
    }


def _panel_4x4(cube: np.ndarray, freqs: np.ndarray, out_png: Path, title_prefix: str) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(13, 12), constrained_layout=True)
    vmin = np.nanpercentile(cube, 2)
    vmax = np.nanpercentile(cube, 98)
    for i, ax in enumerate(axes.ravel()):
        img = cube[:, :, i]
        im = ax.imshow(img, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set_title(f"{freqs[i]:.1f} GHz", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(title_prefix, fontsize=13)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8)
    cbar.ax.tick_params(labelsize=8)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def _decode_if_bytes(v):
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace")
    return str(v)


def _freq_from_id(map_id: str, fallback: float) -> float:
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*GHz", map_id, flags=re.IGNORECASE)
    if m:
        return float(m.group(1))
    return fallback


def _extract_entries_from_map_container(container) -> list[dict]:
    rec = container[0]
    pointer = rec["OMAP"][0]["POINTER"][0]
    ptrs = np.atleast_1d(pointer["PTRS"])
    entries = []
    for i, mp in enumerate(ptrs):
        # IDL list internals may include empty pointer slots; scipy reads these as None.
        if mp is None:
            continue
        map_id = _decode_if_bytes(np.atleast_1d(mp["ID"])[0])
        data_cell = np.asarray(mp["DATA"], dtype=object)
        if data_cell.size == 0:
            continue
        data = np.asarray(data_cell[0], dtype=np.float64)
        freq = _freq_from_id(map_id, float(i))
        stokes = None
        up = map_id.upper()
        if "_I" in up or " I " in up or up.endswith(" I"):
            stokes = "I"
        elif "_V" in up or " V " in up or up.endswith(" V"):
            stokes = "V"
        entries.append({"id": map_id, "data": data, "freq": freq, "stokes": stokes})
    return entries


def _cube_from_entries(entries: list[dict], stokes: str) -> tuple[np.ndarray, np.ndarray]:
    selected = [e for e in entries if e["stokes"] == stokes]
    selected.sort(key=lambda e: e["freq"])
    if not selected:
        raise ValueError(f"No {stokes} maps found in IDL map container.")
    freqs = np.array([e["freq"] for e in selected], dtype=np.float64)
    cube = np.stack([e["data"] for e in selected], axis=-1).astype(np.float64)
    return cube, freqs


def _load_idl_outputs(idl_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    idl = readsav(str(idl_path), verbose=False)

    # Legacy full-state format.
    if "outspace" in idl and "freqlist" in idl:
        ti = np.asarray(idl["outspace"]["TI"][0], dtype=np.float64).transpose(1, 2, 0)
        tv = np.asarray(idl["outspace"]["TV"][0], dtype=np.float64).transpose(1, 2, 0)
        freqs = np.asarray(idl["freqlist"], dtype=np.float64)
        return ti, tv, freqs

    # New compact single-container format.
    if "map" in idl:
        entries = _extract_entries_from_map_container(idl["map"])
        ti, fi = _cube_from_entries(entries, "I")
        tv, fv = _cube_from_entries(entries, "V")
        if ti.shape != tv.shape:
            raise ValueError(f"IDL map container TI/TV shape mismatch: {ti.shape} vs {tv.shape}")
        if fi.shape != fv.shape or not np.allclose(fi, fv, atol=1e-6, rtol=0):
            raise ValueError("IDL map container TI/TV frequency mismatch.")
        return ti, tv, fi

    # Backward compact format.
    if "mapi" in idl and "mapv" in idl:
        ti_entries = _extract_entries_from_map_container(idl["mapi"])
        tv_entries = _extract_entries_from_map_container(idl["mapv"])
        for e in ti_entries:
            e["stokes"] = "I"
        for e in tv_entries:
            e["stokes"] = "V"
        ti, fi = _cube_from_entries(ti_entries, "I")
        tv, fv = _cube_from_entries(tv_entries, "V")
        if ti.shape != tv.shape:
            raise ValueError(f"IDL mapi/mapv TI/TV shape mismatch: {ti.shape} vs {tv.shape}")
        if fi.shape != fv.shape or not np.allclose(fi, fv, atol=1e-6, rtol=0):
            raise ValueError("IDL mapi/mapv TI/TV frequency mismatch.")
        return ti, tv, fi

    raise ValueError(
        f"Unsupported IDL save format in {idl_path}. Expected one of: "
        f"[outspace+freqlist], [map], or [mapi+mapv]."
    )


def _load_python_h5_outputs(py_h5_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(py_h5_path, "r") as f:
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # (nx, ny, nf, 2)
        freqs = np.asarray(f["maps"]["freqlist_ghz"], dtype=np.float64)
    if cube.ndim != 4 or cube.shape[-1] != 2:
        raise ValueError(f"Unexpected Python H5 cube shape: {cube.shape}; expected (nx, ny, nf, 2).")
    # Convert to (ny, nx, nf) to match IDL parser output.
    ti = np.transpose(cube[:, :, :, 0], (1, 0, 2))
    tv = np.transpose(cube[:, :, :, 1], (1, 0, 2))
    return ti, tv, freqs


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare Python-rendered and IDL-rendered MW maps.")
    ap.add_argument("--python-h5", type=Path, default=DEFAULT_PY_H5)
    ap.add_argument("--idl-sav", type=Path, default=DEFAULT_IDL_SAV)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    ti_py, tv_py, freqs_py = _load_python_h5_outputs(args.python_h5)

    ti_idl, tv_idl, freqs_idl = _load_idl_outputs(args.idl_sav)

    if ti_py.shape != ti_idl.shape:
        raise ValueError(f"TI shape mismatch: python={ti_py.shape}, idl={ti_idl.shape}")
    if tv_py.shape != tv_idl.shape:
        raise ValueError(f"TV shape mismatch: python={tv_py.shape}, idl={tv_idl.shape}")
    if not np.allclose(freqs_py, freqs_idl, atol=1e-6, rtol=0):
        raise ValueError(f"Frequency mismatch: python={freqs_py}, idl={freqs_idl}")

    ti_abs = np.abs(ti_py - ti_idl)
    tv_abs = np.abs(tv_py - tv_idl)
    ti_pct = _percent_diff(ti_py, ti_idl)
    tv_pct = _percent_diff(tv_py, tv_idl)
    ti_sym = _sym_diff_pm1(ti_py, ti_idl)
    tv_sym = _sym_diff_pm1(tv_py, tv_idl)
    ti_sum = _sum_diff(ti_py, ti_idl)
    tv_sum = _sum_diff(tv_py, tv_idl)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    _panel_4x4(ti_py, freqs_py, out_dir / "panel_ti_python_4x4.png", "Python TI")
    _panel_4x4(ti_idl, freqs_py, out_dir / "panel_ti_idl_4x4.png", "IDL TI")
    _panel_4x4(ti_pct, freqs_py, out_dir / "panel_ti_pct_diff_python_vs_idl_4x4.png", "TI %DIFF (Python vs IDL)")

    _panel_4x4(tv_py, freqs_py, out_dir / "panel_tv_python_4x4.png", "Python TV")
    _panel_4x4(tv_idl, freqs_py, out_dir / "panel_tv_idl_4x4.png", "IDL TV")
    _panel_4x4(tv_pct, freqs_py, out_dir / "panel_tv_pct_diff_python_vs_idl_4x4.png", "TV %DIFF (Python vs IDL)")
    _panel_4x4(ti_sym, freqs_py, out_dir / "panel_ti_symdiff_pm1_python_vs_idl_4x4.png", "TI SymDiff [-1,1]")
    _panel_4x4(tv_sym, freqs_py, out_dir / "panel_tv_symdiff_pm1_python_vs_idl_4x4.png", "TV SymDiff [-1,1]")

    per_freq = []
    for i, fghz in enumerate(freqs_py):
        per_freq.append(
            {
                "freq_GHz": float(fghz),
                "TI_abs_diff": _stats(ti_abs[:, :, i]),
                "TI_pct_diff_vs_idl": _stats(ti_pct[:, :, i]),
                "TI_sym_diff_pm1": _stats(ti_sym[:, :, i]),
                "TI_sum_diff_py_plus_idl": _stats(ti_sum[:, :, i]),
                "TV_abs_diff": _stats(tv_abs[:, :, i]),
                "TV_pct_diff_vs_idl": _stats(tv_pct[:, :, i]),
                "TV_sym_diff_pm1": _stats(tv_sym[:, :, i]),
                "TV_sum_diff_py_plus_idl": _stats(tv_sum[:, :, i]),
            }
        )

    report = {
        "inputs": {"python_h5": str(args.python_h5), "idl_sav": str(args.idl_sav)},
        "shapes": {"TI": list(ti_py.shape), "TV": list(tv_py.shape)},
        "freqlist_GHz": freqs_py.tolist(),
        "summary": {
            "TI_python": _stats(ti_py),
            "TI_idl": _stats(ti_idl),
            "TI_abs_diff": _stats(ti_abs),
            "TI_pct_diff_vs_idl": _stats(ti_pct),
            "TI_sym_diff_pm1": _stats(ti_sym),
            "TI_sum_diff_py_plus_idl": _stats(ti_sum),
            "TV_python": _stats(tv_py),
            "TV_idl": _stats(tv_idl),
            "TV_abs_diff": _stats(tv_abs),
            "TV_pct_diff_vs_idl": _stats(tv_pct),
            "TV_sym_diff_pm1": _stats(tv_sym),
            "TV_sum_diff_py_plus_idl": _stats(tv_sum),
        },
        "per_frequency": per_freq,
    }

    with (out_dir / "comparison_python_vs_idl.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("Outputs:")
    print(f"- comparison_json: {out_dir / 'comparison_python_vs_idl.json'}")
    print(f"- panel_ti_pct_diff: {out_dir / 'panel_ti_pct_diff_python_vs_idl_4x4.png'}")
    print(f"- panel_tv_pct_diff: {out_dir / 'panel_tv_pct_diff_python_vs_idl_4x4.png'}")
    print(f"- panel_ti_symdiff: {out_dir / 'panel_ti_symdiff_pm1_python_vs_idl_4x4.png'}")
    print(f"- panel_tv_symdiff: {out_dir / 'panel_tv_symdiff_pm1_python_vs_idl_4x4.png'}")
    print("Summary TI %diff mean/std/min/max:", report["summary"]["TI_pct_diff_vs_idl"])
    print("Summary TV %diff mean/std/min/max:", report["summary"]["TV_pct_diff_vs_idl"])
    print("Summary TI sym(-1..1) mean/std/min/max:", report["summary"]["TI_sym_diff_pm1"])
    print("Summary TV sym(-1..1) mean/std/min/max:", report["summary"]["TV_sym_diff_pm1"])


if __name__ == "__main__":
    main()
