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


def _panel_4x4(
    cube: np.ndarray,
    plane_vals: np.ndarray,
    out_png: Path,
    title_prefix: str,
    *,
    plane_unit: str = "GHz",
) -> None:
    fig, axes = plt.subplots(4, 4, figsize=(13, 12), constrained_layout=True)
    finite_cube = np.asarray(cube, dtype=np.float64)
    finite_vals = finite_cube[np.isfinite(finite_cube)]
    if finite_vals.size:
        vmin = float(np.nanpercentile(finite_vals, 2))
        vmax = float(np.nanpercentile(finite_vals, 98))
    else:
        vmin, vmax = 0.0, 1.0
    im = None
    nplanes = int(cube.shape[2]) if cube.ndim == 3 else 0
    for i, ax in enumerate(axes.ravel()):
        if i >= nplanes:
            ax.axis("off")
            continue
        img = cube[:, :, i]
        im = ax.imshow(img, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        if plane_unit:
            ax.set_title(f"{plane_vals[i]:.1f} {plane_unit}", fontsize=9)
        else:
            v = plane_vals[i]
            if float(v).is_integer():
                ax.set_title(f"{int(v)}", fontsize=9)
            else:
                ax.set_title(f"{v:.3g}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(title_prefix, fontsize=13)
    if im is not None:
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
    # EUV labels often look like "... A94" or "... 171"
    m = re.search(r"\bA?([0-9]{2,4})\b", map_id, flags=re.IGNORECASE)
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
        # Reuse the MW comparison pipeline for EUV by mapping:
        # CORONA -> I, TR -> V
        elif "GX (CORONA)" in up:
            stokes = "I"
        elif "GX (TR)" in up:
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


def _load_python_h5_outputs(py_h5_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    with h5py.File(py_h5_path, "r") as f:
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # (nx, ny, nf, 2)
        schema = "mw" if "freqlist_ghz" in f["maps"] else ("euv" if "channel_ids" in f["maps"] else "unknown")
        if "freqlist_ghz" in f["maps"]:
            freqs = np.asarray(f["maps"]["freqlist_ghz"], dtype=np.float64)
            meta = {
                "kind": "mw",
                "comp1_label": "TI",
                "comp2_label": "TV",
                "plane_unit": "GHz",
                "plane_name": "freq",
                "plane_ids": [float(x) for x in np.asarray(freqs).tolist()],
            }
        elif "channel_ids" in f["maps"]:
            raw = f["maps"]["channel_ids"][()]
            comp_raw = f["maps"]["component_ids"][()] if "component_ids" in f["maps"] else np.asarray([b"COMP1", b"COMP2"])
            vals = []
            ids = []
            for x in np.asarray(raw).reshape(-1):
                s = x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)
                ids.append(s)
                m = re.search(r"([0-9]+(?:\.[0-9]+)?)", s)
                vals.append(float(m.group(1)) if m else float(len(vals)))
            freqs = np.asarray(vals, dtype=np.float64)
            comp_ids = [
                (x.decode("utf-8", errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)).upper()
                for x in np.asarray(comp_raw).reshape(-1)
            ]
            comp1 = comp_ids[0] if len(comp_ids) > 0 else "COMP1"
            comp2 = comp_ids[1] if len(comp_ids) > 1 else "COMP2"
            meta = {
                "kind": "euv",
                "comp1_label": comp1,
                "comp2_label": comp2,
                "plane_unit": "",
                "plane_name": "channel",
                "plane_ids": ids,
            }
        else:
            raise ValueError("Unsupported Python H5 schema: expected maps/freqlist_ghz or maps/channel_ids.")
    if cube.ndim != 4 or cube.shape[-1] != 2:
        raise ValueError(f"Unexpected Python H5 cube shape: {cube.shape}; expected (nx, ny, nf, 2).")
    # Convert to (ny, nx, nf) to match IDL parser output.
    ti = np.transpose(cube[:, :, :, 0], (1, 0, 2))
    tv = np.transpose(cube[:, :, :, 1], (1, 0, 2))
    if schema not in {"mw", "euv"}:
        meta = {
            "kind": "mw",
            "comp1_label": "TI",
            "comp2_label": "TV",
            "plane_unit": "GHz",
            "plane_name": "freq",
            "plane_ids": [float(x) for x in np.asarray(freqs).tolist()],
        }
    return ti, tv, freqs, meta


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare Python-rendered and IDL-rendered MW/EUV maps.")
    ap.add_argument("--python-h5", type=Path, default=DEFAULT_PY_H5)
    ap.add_argument("--idl-sav", type=Path, default=DEFAULT_IDL_SAV)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument(
        "--kind",
        choices=["auto", "mw", "euv"],
        default="auto",
        help="Comparison labeling mode. 'auto' infers from the Python H5 schema.",
    )
    args = ap.parse_args()

    ti_py, tv_py, freqs_py, py_meta = _load_python_h5_outputs(args.python_h5)
    kind = py_meta["kind"] if args.kind == "auto" else args.kind
    comp1_label = "TI" if kind == "mw" else "CORONA"
    comp2_label = "TV" if kind == "mw" else "TR"
    plane_unit = "GHz" if kind == "mw" else ""
    if kind == "euv" and py_meta.get("plane_ids"):
        plane_vals_for_plot = freqs_py
    else:
        plane_vals_for_plot = freqs_py

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

    comp1_slug = comp1_label.lower()
    comp2_slug = comp2_label.lower()
    panel_comp1_pct = out_dir / f"panel_{comp1_slug}_pct_diff_python_vs_idl_4x4.png"
    panel_comp2_pct = out_dir / f"panel_{comp2_slug}_pct_diff_python_vs_idl_4x4.png"
    panel_comp1_sym = out_dir / f"panel_{comp1_slug}_symdiff_pm1_python_vs_idl_4x4.png"
    panel_comp2_sym = out_dir / f"panel_{comp2_slug}_symdiff_pm1_python_vs_idl_4x4.png"

    _panel_4x4(ti_py, plane_vals_for_plot, out_dir / f"panel_{comp1_slug}_python_4x4.png", f"Python {comp1_label}", plane_unit=plane_unit)
    _panel_4x4(ti_idl, plane_vals_for_plot, out_dir / f"panel_{comp1_slug}_idl_4x4.png", f"IDL {comp1_label}", plane_unit=plane_unit)
    _panel_4x4(ti_pct, plane_vals_for_plot, panel_comp1_pct, f"{comp1_label} %DIFF (Python vs IDL)", plane_unit=plane_unit)

    _panel_4x4(tv_py, plane_vals_for_plot, out_dir / f"panel_{comp2_slug}_python_4x4.png", f"Python {comp2_label}", plane_unit=plane_unit)
    _panel_4x4(tv_idl, plane_vals_for_plot, out_dir / f"panel_{comp2_slug}_idl_4x4.png", f"IDL {comp2_label}", plane_unit=plane_unit)
    _panel_4x4(tv_pct, plane_vals_for_plot, panel_comp2_pct, f"{comp2_label} %DIFF (Python vs IDL)", plane_unit=plane_unit)
    _panel_4x4(ti_sym, plane_vals_for_plot, panel_comp1_sym, f"{comp1_label} SymDiff [-1,1]", plane_unit=plane_unit)
    _panel_4x4(tv_sym, plane_vals_for_plot, panel_comp2_sym, f"{comp2_label} SymDiff [-1,1]", plane_unit=plane_unit)

    per_freq = []
    for i, fghz in enumerate(freqs_py):
        per_freq.append(
            {
                ("freq_GHz" if kind == "mw" else "channel"): (float(fghz) if kind == "mw" else py_meta["plane_ids"][i]),
                f"{comp1_label}_abs_diff": _stats(ti_abs[:, :, i]),
                f"{comp1_label}_pct_diff_vs_idl": _stats(ti_pct[:, :, i]),
                f"{comp1_label}_sym_diff_pm1": _stats(ti_sym[:, :, i]),
                f"{comp1_label}_sum_diff_py_plus_idl": _stats(ti_sum[:, :, i]),
                f"{comp2_label}_abs_diff": _stats(tv_abs[:, :, i]),
                f"{comp2_label}_pct_diff_vs_idl": _stats(tv_pct[:, :, i]),
                f"{comp2_label}_sym_diff_pm1": _stats(tv_sym[:, :, i]),
                f"{comp2_label}_sum_diff_py_plus_idl": _stats(tv_sum[:, :, i]),
            }
        )

    report = {
        "inputs": {"python_h5": str(args.python_h5), "idl_sav": str(args.idl_sav)},
        "mode": {
            "kind": kind,
            "comp1_label": comp1_label,
            "comp2_label": comp2_label,
            "plane_unit": plane_unit,
            "plane_ids": py_meta.get("plane_ids", [float(x) for x in freqs_py.tolist()]),
        },
        "shapes": {comp1_label: list(ti_py.shape), comp2_label: list(tv_py.shape)},
        ("freqlist_GHz" if kind == "mw" else "channels"): (freqs_py.tolist() if kind == "mw" else py_meta.get("plane_ids", [])),
        "summary": {
            f"{comp1_label}_python": _stats(ti_py),
            f"{comp1_label}_idl": _stats(ti_idl),
            f"{comp1_label}_abs_diff": _stats(ti_abs),
            f"{comp1_label}_pct_diff_vs_idl": _stats(ti_pct),
            f"{comp1_label}_sym_diff_pm1": _stats(ti_sym),
            f"{comp1_label}_sum_diff_py_plus_idl": _stats(ti_sum),
            f"{comp2_label}_python": _stats(tv_py),
            f"{comp2_label}_idl": _stats(tv_idl),
            f"{comp2_label}_abs_diff": _stats(tv_abs),
            f"{comp2_label}_pct_diff_vs_idl": _stats(tv_pct),
            f"{comp2_label}_sym_diff_pm1": _stats(tv_sym),
            f"{comp2_label}_sum_diff_py_plus_idl": _stats(tv_sum),
        },
        "per_frequency": per_freq,
    }

    with (out_dir / "comparison_python_vs_idl.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("Outputs:")
    print(f"- comparison_json: {out_dir / 'comparison_python_vs_idl.json'}")
    print(f"- panel_{comp1_slug}_pct_diff: {panel_comp1_pct}")
    print(f"- panel_{comp2_slug}_pct_diff: {panel_comp2_pct}")
    print(f"- panel_{comp1_slug}_symdiff: {panel_comp1_sym}")
    print(f"- panel_{comp2_slug}_symdiff: {panel_comp2_sym}")
    print(f"Summary {comp1_label} %diff mean/std/min/max:", report["summary"][f"{comp1_label}_pct_diff_vs_idl"])
    print(f"Summary {comp2_label} %diff mean/std/min/max:", report["summary"][f"{comp2_label}_pct_diff_vs_idl"])
    print(f"Summary {comp1_label} sym(-1..1) mean/std/min/max:", report["summary"][f"{comp1_label}_sym_diff_pm1"])
    print(f"Summary {comp2_label} sym(-1..1) mean/std/min/max:", report["summary"][f"{comp2_label}_sym_diff_pm1"])


if __name__ == "__main__":
    main()
