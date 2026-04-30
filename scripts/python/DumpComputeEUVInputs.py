#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

# Allow this test utility to run in restricted environments without requiring
# callers to preconfigure SunPy writable paths.
if not os.environ.get("SUNPY_CONFIGDIR"):
    _sunpy_cfg = Path("/tmp/.sunpy-config")
    _sunpy_cfg.mkdir(parents=True, exist_ok=True)
    os.environ["SUNPY_CONFIGDIR"] = str(_sunpy_cfg)

from gxrender.euv import GXEUVImageComputing, load_euv_response_sav
from gxrender.utils.test_data import find_default_model_file, test_data_setup_hint, try_find_response_file
from gxrender.workflows._render_common import DEFAULT_OUTDIR, prepare_common_inputs


def _example_default_shtable() -> np.ndarray:
    weights = np.array([1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4], dtype=np.float64)
    shtable = np.outer(weights, weights)
    shtable[6, 6] = 0.1
    return shtable


def _example_plasma_defaults() -> dict[str, float]:
    return {
        "tbase": 1e6,
        "nbase": 1e8,
        "q0": 0.0217,
        "a": 0.3,
        "b": 2.7,
    }


def _build_coronaparms_dtype() -> np.dtype:
    return np.dtype(
        [
            ("Tbase", np.float64),
            ("nbase", np.float64),
            ("Q0", np.float64),
            ("a", np.float64),
            ("b", np.float64),
            ("mode", np.int32),
        ]
    )


def _resolve_default_response_sav(instrument: str) -> Path | None:
    env_path = os.environ.get("GXIMAGECOMPUTING_EUV_RESPONSE_SAV", "").strip()
    if env_path:
        p = Path(env_path).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"GXIMAGECOMPUTING_EUV_RESPONSE_SAV points to a missing file: {p}")

    repo_response = try_find_response_file(instrument)
    if repo_response is not None:
        return repo_response
    return None


def _build_coronaparms(
    *,
    tbase: float,
    nbase: float,
    q0: float,
    a: float,
    b: float,
    mode: int,
) -> tuple[np.ndarray, np.dtype]:
    dt_c = _build_coronaparms_dtype()
    cparms = np.zeros(1, dtype=dt_c)
    cparms["Tbase"] = float(tbase)
    cparms["nbase"] = float(nbase)
    cparms["Q0"] = float(q0)
    cparms["a"] = float(a)
    cparms["b"] = float(b)
    cparms["mode"] = int(mode)
    return cparms, dt_c


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Dump exact ComputeEUV inputs (Python side) before the DLL call.")
    p.add_argument("--model-path", type=Path, default=None)
    p.add_argument("--model-format", choices=["h5", "sav", "auto"], default="auto")
    p.add_argument(
        "--ebtel-path",
        type=str,
        default=None,
        help='Optional EBTEL table (.sav). Use "" to disable.',
    )
    p.add_argument(
        "--response-sav",
        type=Path,
        default=None,
        help="Path to IDL-like EUV response SAV.",
    )
    p.add_argument(
        "--channels",
        nargs="*",
        default=["94", "131", "171", "193", "211", "304", "335"],
    )
    p.add_argument("--instrument", type=str, default="AIA")
    p.add_argument("--omp-threads", type=int, default=8)
    p.add_argument("--xc", type=float, default=None)
    p.add_argument("--yc", type=float, default=None)
    p.add_argument(
        "--dsun-cm",
        type=float,
        default=None,
        help="Override model.DSun before preparing DLL inputs (cm).",
    )
    p.add_argument(
        "--lonc-deg",
        type=float,
        default=None,
        help="Override model.lonC before preparing DLL inputs (deg).",
    )
    p.add_argument(
        "--b0sun-deg",
        type=float,
        default=None,
        help="Override model.b0Sun before preparing DLL inputs (deg).",
    )
    p.add_argument("--dx", type=float, default=None)
    p.add_argument("--dy", type=float, default=None)
    p.add_argument("--pixel-scale-arcsec", type=float, default=None)
    p.add_argument("--nx", type=int, default=None)
    p.add_argument("--ny", type=int, default=None)
    p.add_argument(
        "--xrange", type=float, nargs=2, default=None, metavar=("XMIN", "XMAX")
    )
    p.add_argument(
        "--yrange", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX")
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUTDIR / "computeeuv_inputs_python",
        help="Directory for .npy dumps + manifest.json",
    )
    p.add_argument(
        "--mode",
        type=int,
        default=0,
        help="Corona params mode flag (passed to ComputeEUV).",
    )
    return p.parse_args()


def _save_array(out_dir: Path, name: str, arr: np.ndarray) -> None:
    np.save(out_dir / f"{name}.npy", arr, allow_pickle=False)


def main() -> None:
    args = _parse_args()
    if args.model_path is None:
        args.model_path = find_default_model_file(".h5")
    if args.ebtel_path is None:
        args.ebtel_path = ""
    if args.pixel_scale_arcsec is None and args.dx is None and args.dy is None:
        args.pixel_scale_arcsec = 2.0
    common = prepare_common_inputs(args, prefer_execute_center=False)

    if args.response_sav is not None:
        response_path = str(args.response_sav)
        response, response_dt, response_meta = load_euv_response_sav(response_path)
        response_source_mode = "explicit"
    else:
        auto_response_sav = _resolve_default_response_sav(args.instrument)
        if auto_response_sav is not None:
            response_path = str(auto_response_sav)
            response, response_dt, response_meta = load_euv_response_sav(response_path)
            response_source_mode = "auto_sav"
        else:
            raise FileNotFoundError(
                "No explicit EUV response SAV was provided, and no default response fixture could be found. "
                + test_data_setup_hint(f"EUV response file for instrument {str(args.instrument).strip().lower()!r}")
                + " You may also set GXIMAGECOMPUTING_EUV_RESPONSE_SAV to an explicit SAV file."
            )

    plasma = _example_plasma_defaults()
    coronaparms, coronaparms_dt = _build_coronaparms(
        tbase=plasma["tbase"],
        nbase=plasma["nbase"],
        q0=plasma["q0"],
        a=plasma["a"],
        b=plasma["b"],
        mode=int(args.mode),
    )

    gx = GXEUVImageComputing()
    simbox, simbox_dt = gx._make_euv_simbox(
        box_nx=common.nx,
        box_ny=common.ny,
        box_xc=common.xc,
        box_yc=common.yc,
        box_dx=common.dx,
        box_dy=common.dy,
        parallel=False,
        exact=False,
        nthreads=0,
    )
    outspace, outspace_dt = gx._reserve_euv_output(
        box_nx=common.nx,
        box_ny=common.ny,
        nchan=int(response["Nchannels"][0]),
    )

    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    _save_array(out_dir, "model", np.asarray(common.model))
    _save_array(out_dir, "ebtel", np.asarray(common.ebtel_c))
    _save_array(out_dir, "response", np.asarray(response))
    _save_array(out_dir, "simbox", np.asarray(simbox))
    _save_array(out_dir, "coronaparms", np.asarray(coronaparms))
    _save_array(out_dir, "outspace", np.asarray(outspace))
    _save_array(out_dir, "shtable", _example_default_shtable())

    manifest = {
        "inputs": {
            "model_path": str(common.model_path),
            "model_format": str(common.loader),
            "ebtel_path": str(common.ebtel_path),
            "response_sav": response_path,
            "response_source_mode": response_source_mode,
            "instrument": str(response_meta.instrument),
            "channels": [str(c) for c in response_meta.channels],
        },
        "geometry": {
            "center_source": str(common.center_source),
            "xc_arcsec": float(common.xc),
            "yc_arcsec": float(common.yc),
            "dx_arcsec": float(common.dx),
            "dy_arcsec": float(common.dy),
            "nx": int(common.nx),
            "ny": int(common.ny),
            "fov_x_arcsec": float(common.fov_x),
            "fov_y_arcsec": float(common.fov_y),
        },
        "observer_overrides_applied": common.observer_overrides_applied,
        "dll_inputs": {
            "model_dtype": np.asarray(common.model).dtype.descr,
            "ebtel_dtype": np.asarray(common.ebtel_c).dtype.descr,
            "response_dtype": np.asarray(response).dtype.descr,
            "simbox_dtype": np.asarray(simbox).dtype.descr,
            "coronaparms_dtype": np.asarray(coronaparms).dtype.descr,
            "outspace_dtype": np.asarray(outspace).dtype.descr,
            "shtable_shape": [7, 7],
        },
        "internal_dtypes": {
            "response_dt": response_dt.descr,
            "simbox_dt": simbox_dt.descr,
            "coronaparms_dt": coronaparms_dt.descr,
            "outspace_dt": outspace_dt.descr,
        },
        "notes": "Arrays are dumped exactly before the Python ComputeEUV DLL call.",
    }
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as fp:
        json.dump(manifest, fp, indent=2)

    print("Dumped Python ComputeEUV inputs:")
    print(f"- out_dir: {out_dir}")
    for name in (
        "model",
        "ebtel",
        "response",
        "simbox",
        "coronaparms",
        "outspace",
        "shtable",
    ):
        arr = np.load(out_dir / f"{name}.npy", allow_pickle=False)
        print(f"- {name}: shape={arr.shape} dtype={arr.dtype}")
    print(f"- manifest: {out_dir / 'manifest.json'}")


if __name__ == "__main__":
    main()
