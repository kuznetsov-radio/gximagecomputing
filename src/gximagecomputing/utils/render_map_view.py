from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.widgets import CheckButtons, RangeSlider, Slider


def _decode_h5_scalar(v) -> str:
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace")
    return str(v)


def _read_render_h5(path: Path):
    with h5py.File(path, "r") as f:
        if "maps" not in f or "data" not in f["maps"]:
            raise ValueError(f"{path} is not a rendered map H5 (missing maps/data)")
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # [nx, ny, nf, 2]
        freqs = np.asarray(f["maps"]["freqlist_ghz"], dtype=np.float64)
        if cube.ndim != 4 or cube.shape[-1] != 2:
            raise ValueError(f"Expected cube shape [nx, ny, nf, 2], got {cube.shape}")
        if cube.shape[2] != freqs.shape[0]:
            raise ValueError(
                f"Frequency axis mismatch: cube nf={cube.shape[2]} vs freqlist={freqs.shape[0]}"
            )

        index_header = ""
        if "metadata" in f and "index_header" in f["metadata"]:
            index_header = _decode_h5_scalar(f["metadata"]["index_header"][()])

        date_obs = ""
        if "metadata" in f and "date_obs" in f["metadata"]:
            date_obs = _decode_h5_scalar(f["metadata"]["date_obs"][()])

    return cube, freqs, index_header, date_obs


def _norm_for_data(data: np.ndarray, vmin: float, vmax: float, log: bool):
    if log:
        vvmin = max(vmin, np.finfo(float).tiny)
        vvmax = max(vmax, vvmin * 1.000001)
        return LogNorm(vmin=vvmin, vmax=vvmax)
    if vmin < 0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    return Normalize(vmin=vmin, vmax=vmax)


def _safe_minmax(a: np.ndarray) -> tuple[float, float]:
    finite = np.isfinite(a)
    if not finite.any():
        return 0.0, 1.0
    mn = float(np.nanmin(a))
    mx = float(np.nanmax(a))
    if not np.isfinite(mn) or not np.isfinite(mx):
        return 0.0, 1.0
    if mx <= mn:
        mx = mn + 1.0
    return mn, mx


def _title(path: Path, stokes: str, freq: float, date_obs: str) -> str:
    dt = f"\n{date_obs}" if date_obs else ""
    return f"{path.name} | {stokes} @ {freq:.2f} GHz{dt}"


def run_viewer(path: Path, start_index: int = 0) -> None:
    cube, freqs, _, date_obs = _read_render_h5(path)

    # Convert to matplotlib-friendly display axes: [ny, nx, nf, 2]
    disp = np.transpose(cube, (1, 0, 2, 3))
    ti_cube = disp[:, :, :, 0]
    tv_cube = disp[:, :, :, 1]
    nf = ti_cube.shape[2]

    idx0 = int(np.clip(start_index, 0, nf - 1))

    ti_min, ti_max = _safe_minmax(ti_cube)
    tv_min, tv_max = _safe_minmax(tv_cube)

    fig = plt.figure(figsize=(12, 7.8))
    # Manual layout for stable spacing across backends.
    # Top control bar
    cb_ax = fig.add_axes([0.06, 0.93, 0.10, 0.05], frameon=False)  # TI log
    cb_ax.set_xticks([])
    cb_ax.set_yticks([])
    ax_freq = fig.add_axes([0.24, 0.94, 0.62, 0.03])  # frequency slider
    # Main map panels
    ax_ti = fig.add_axes([0.06, 0.27, 0.38, 0.60])
    ax_tv = fig.add_axes([0.56, 0.27, 0.38, 0.60])
    # Range sliders: each constrained to its map width, separated in the middle.
    ax_ti_rng = fig.add_axes([0.06, 0.14, 0.36, 0.035])
    ax_tv_rng = fig.add_axes([0.58, 0.14, 0.36, 0.035])
    # Dedicated range-info labels to avoid overlap with RangeSlider default valtext.
    ti_info = fig.text(0.06, 0.182, "", ha="left", va="bottom")
    tv_info = fig.text(0.58, 0.182, "", ha="left", va="bottom")

    ti_data = ti_cube[:, :, idx0]
    tv_data = tv_cube[:, :, idx0]

    ti_norm = _norm_for_data(ti_data, ti_min, ti_max, log=False)
    tv_norm = _norm_for_data(tv_data, tv_min, tv_max, log=False)

    im_ti = ax_ti.imshow(ti_data, origin="lower", cmap="inferno", norm=ti_norm)
    im_tv = ax_tv.imshow(tv_data, origin="lower", cmap="coolwarm", norm=tv_norm)

    ax_ti.set_title(_title(path, "TI", freqs[idx0], date_obs))
    ax_tv.set_title(_title(path, "TV", freqs[idx0], date_obs))
    ax_ti.set_xlabel("X [pix]")
    ax_ti.set_ylabel("Y [pix]")
    ax_tv.set_xlabel("X [pix]")
    ax_tv.set_ylabel("Y [pix]")

    cbar_ti = fig.colorbar(im_ti, ax=ax_ti)
    cbar_tv = fig.colorbar(im_tv, ax=ax_tv)
    cbar_ti.set_label("Brightness Temperature [K]")
    cbar_tv.set_label("Brightness Temperature [K]")

    freq_slider = Slider(ax=ax_freq, label="Frequency index", valmin=0, valmax=nf - 1, valinit=idx0, valstep=1)
    ti_range = RangeSlider(
        ax=ax_ti_rng,
        label="",
        valmin=ti_min,
        valmax=ti_max,
        valinit=(ti_min, ti_max),
    )
    tv_range = RangeSlider(
        ax=ax_tv_rng,
        label="",
        valmin=tv_min,
        valmax=tv_max,
        valinit=(tv_min, tv_max),
    )
    ti_range.valtext.set_visible(False)
    tv_range.valtext.set_visible(False)

    # TI log in top control bar (left), aligned with frequency slider.
    cb = CheckButtons(cb_ax, ["TI log"], [False])

    state = {"idx": idx0, "ti_log": False}

    def redraw() -> None:
        idx = state["idx"]
        ti = ti_cube[:, :, idx]
        tv = tv_cube[:, :, idx]

        im_ti.set_data(ti)
        im_tv.set_data(tv)

        ti_vmin, ti_vmax = ti_range.val
        tv_vmin, tv_vmax = tv_range.val

        im_ti.set_norm(_norm_for_data(ti, float(ti_vmin), float(ti_vmax), state["ti_log"]))
        im_tv.set_norm(_norm_for_data(tv, float(tv_vmin), float(tv_vmax), False))
        ti_info.set_text(f"TI range [K]: ({ti_vmin:.3g}, {ti_vmax:.3g})")
        tv_info.set_text(f"TV range [K]: ({tv_vmin:.3g}, {tv_vmax:.3g})")

        ax_ti.set_title(_title(path, "TI", freqs[idx], date_obs))
        ax_tv.set_title(_title(path, "TV", freqs[idx], date_obs))

        fig.canvas.draw_idle()

    def on_freq(val):
        state["idx"] = int(val)
        redraw()

    def on_ti_range(_):
        redraw()

    def on_tv_range(_):
        redraw()

    def on_log(_):
        state["ti_log"] = not state["ti_log"]
        redraw()

    freq_slider.on_changed(on_freq)
    ti_range.on_changed(on_ti_range)
    tv_range.on_changed(on_tv_range)
    cb.on_clicked(on_log)
    redraw()

    plt.show()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Interactive viewer for RenderExampleMW.py output H5 files ([nx, ny, nf, 2])."
    )
    p.add_argument("h5_path", type=Path, help="Path to rendered map H5.")
    p.add_argument("--start-index", type=int, default=0, help="Initial frequency index.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_viewer(args.h5_path.expanduser().resolve(), start_index=args.start_index)


if __name__ == "__main__":
    main()
