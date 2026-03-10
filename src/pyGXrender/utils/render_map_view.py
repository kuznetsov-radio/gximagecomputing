from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.io import fits
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.widgets import CheckButtons, RangeSlider, Slider
from scipy.io import readsav


class _ViewerData(dict):
    pass


def _decode_scalar(v):
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace")
    return str(v)


def _extract_freq_from_id(map_id: str, fallback: float) -> float:
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*GHz", map_id, flags=re.IGNORECASE)
    return float(m.group(1)) if m else float(fallback)


def _extract_channel_from_id(map_id: str, fallback: str) -> str:
    m = re.search(r"\b(?:AIA\s+)?A?(94|131|171|193|211|304|335)\b", map_id, flags=re.IGNORECASE)
    return m.group(1) if m else fallback


def _extract_entries_from_map_container(container) -> list[dict]:
    rec = container[0]
    pointer = rec["OMAP"][0]["POINTER"][0]
    ptrs = np.atleast_1d(pointer["PTRS"])
    entries: list[dict] = []
    for i, mp in enumerate(ptrs):
        if mp is None:
            continue
        data_cell = np.asarray(mp["DATA"], dtype=object)
        if data_cell.size == 0:
            continue
        data = np.asarray(data_cell[0], dtype=np.float64)
        map_id = _decode_scalar(np.atleast_1d(mp["ID"])[0])
        entries.append({"id": map_id, "data": data, "index": int(i)})
    return entries


def _classify_combined_entries(entries: list[dict]) -> tuple[list[dict], list[dict], str | None]:
    mw_i, mw_v, euv_tr, euv_cor = [], [], [], []
    for e in entries:
        up = e["id"].upper()
        e["freq"] = _extract_freq_from_id(e["id"], e["index"])
        e["channel"] = _extract_channel_from_id(e["id"], str(e["index"]))
        if "TB_I" in up or "_I" in up or up.endswith(" I") or " I " in up:
            mw_i.append(e)
        elif "TB_V" in up or "_V" in up or up.endswith(" V") or " V " in up:
            mw_v.append(e)
        elif "(TR)" in up or " TR" in up:
            euv_tr.append(e)
        elif "(CORONA)" in up or "CORONA" in up:
            euv_cor.append(e)
    if mw_i and mw_v:
        return mw_i, mw_v, "mw"
    if euv_tr and euv_cor:
        return euv_tr, euv_cor, "euv"
    return [], [], None


def _stack_entries_by_freq(entries: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    items = sorted(entries, key=lambda e: e["freq"])
    return (
        np.stack([np.asarray(e["data"], dtype=np.float64) for e in items], axis=-1),
        np.asarray([float(e["freq"]) for e in items], dtype=np.float64),
    )


def _stack_entries_by_channel(entries: list[dict], channel_order: list[str] | None = None) -> tuple[np.ndarray, list[str]]:
    grouped = {str(e["channel"]): e for e in entries}
    if channel_order is None:
        keys = list(grouped.keys())
        try:
            channel_order = sorted(keys, key=lambda s: float(s))
        except Exception:
            channel_order = keys
    cube = np.stack([np.asarray(grouped[ch]["data"], dtype=np.float64) for ch in channel_order], axis=-1)
    return cube, list(channel_order)


def _read_render_h5(path: Path) -> _ViewerData:
    with h5py.File(path, "r") as f:
        if "maps" not in f or "data" not in f["maps"]:
            raise ValueError(f"{path} is not a rendered map H5 (missing maps/data)")
        cube = np.asarray(f["maps"]["data"], dtype=np.float64)  # [nx, ny, n, 2]
        if cube.ndim != 4 or cube.shape[-1] != 2:
            raise ValueError(f"Expected cube shape [nx, ny, n, 2], got {cube.shape}")

        maps = f["maps"]
        meta = f.get("metadata")
        index_header = _decode_scalar(meta["index_header"][()]) if meta is not None and "index_header" in meta else ""
        date_obs = _decode_scalar(meta["date_obs"][()]) if meta is not None and "date_obs" in meta else ""

        bunit = ""
        if index_header:
            try:
                hdr = fits.Header.fromstring(index_header, sep="\n")
                bunit = str(hdr.get("BUNIT", ""))
            except Exception:
                pass

        disp = np.transpose(cube, (1, 0, 2, 3))
        left_cube = disp[:, :, :, 0]
        right_cube = disp[:, :, :, 1]

        if "freqlist_ghz" in maps:
            freqs = np.asarray(maps["freqlist_ghz"], dtype=np.float64)
            if left_cube.shape[2] != freqs.shape[0]:
                raise ValueError(f"Frequency axis mismatch: cube n={left_cube.shape[2]} vs freqlist={freqs.shape[0]}")
            return _ViewerData(
                kind="mw",
                left_cube=left_cube,
                right_cube=right_cube,
                axis_kind="freq",
                axis_values=freqs,
                axis_labels=[f"{f:.2f}" for f in freqs],
                left_label="TI",
                right_label="TV",
                left_cmap="inferno",
                right_cmap="coolwarm",
                bunit=bunit or "K",
                date_obs=date_obs,
                index_header=index_header,
            )

        if "channel_ids" in maps:
            channels = [_decode_scalar(v) for v in np.asarray(maps["channel_ids"])]
            comp = ["TR", "Corona"]
            if "component_ids" in maps:
                comp = [_decode_scalar(v) for v in np.asarray(maps["component_ids"])[:2]]
            if left_cube.shape[2] != len(channels):
                raise ValueError(f"Channel axis mismatch: cube n={left_cube.shape[2]} vs channels={len(channels)}")
            comp_up = [c.upper() for c in comp]
            # Normalize EUV display order to TR (left), Corona (right), even if stored as [Corona, TR].
            if len(comp_up) >= 2 and comp_up[0].startswith("COR") and comp_up[1].startswith("TR"):
                left_cube, right_cube = right_cube, left_cube
                comp = [comp[1], comp[0]]
                comp_up = [c.upper() for c in comp]
            return _ViewerData(
                kind="euv",
                left_cube=left_cube,
                right_cube=right_cube,
                axis_kind="channel",
                axis_values=np.arange(len(channels), dtype=np.float64),
                axis_labels=channels,
                left_label=f"GX ({'TR' if comp_up[0].startswith('TR') else 'Corona'})",
                right_label=f"GX ({'TR' if comp_up[1].startswith('TR') else 'Corona'})",
                left_cmap="magma",
                right_cmap="magma",
                bunit=bunit or "DN s^-1 pix^-1",
                date_obs=date_obs,
                index_header=index_header,
            )

    raise ValueError(f"{path} H5 maps container is neither MW (freqlist_ghz) nor EUV (channel_ids).")


def _read_render_sav(path: Path) -> _ViewerData:
    idl = readsav(str(path), verbose=False)
    key_lookup = {k.lower(): k for k in idl.keys()}

    if "map" in key_lookup:
        entries = _extract_entries_from_map_container(idl[key_lookup["map"]])
        left_entries, right_entries, kind = _classify_combined_entries(entries)
        if kind == "mw":
            left_cube, freqs = _stack_entries_by_freq(left_entries)
            right_cube, freqs2 = _stack_entries_by_freq(right_entries)
            if left_cube.shape != right_cube.shape:
                raise ValueError(f"MW SAV TI/TV shape mismatch: {left_cube.shape} vs {right_cube.shape}")
            if freqs.shape != freqs2.shape or not np.allclose(freqs, freqs2, atol=1e-6, rtol=0):
                raise ValueError("MW SAV TI/TV frequency mismatch")
            return _ViewerData(
                kind="mw",
                left_cube=left_cube,
                right_cube=right_cube,
                axis_kind="freq",
                axis_values=freqs,
                axis_labels=[f"{f:.2f}" for f in freqs],
                left_label="TI",
                right_label="TV",
                left_cmap="inferno",
                right_cmap="coolwarm",
                bunit="K",
                date_obs="",
                index_header="",
            )
        if kind == "euv":
            tr_cube, channels = _stack_entries_by_channel(left_entries)
            cor_cube, channels2 = _stack_entries_by_channel(right_entries, channel_order=channels)
            if tr_cube.shape != cor_cube.shape:
                raise ValueError(f"EUV SAV TR/Corona shape mismatch: {tr_cube.shape} vs {cor_cube.shape}")
            return _ViewerData(
                kind="euv",
                left_cube=tr_cube,
                right_cube=cor_cube,
                axis_kind="channel",
                axis_values=np.arange(len(channels), dtype=np.float64),
                axis_labels=channels2,
                left_label="GX (TR)",
                right_label="GX (Corona)",
                left_cmap="magma",
                right_cmap="magma",
                bunit="DN s^-1 pix^-1",
                date_obs="",
                index_header="",
            )
        raise ValueError(f"Unsupported combined IDL map container in {path}")

    if "mapcorona" in key_lookup and "maptr" in key_lookup:
        cor_entries = _extract_entries_from_map_container(idl[key_lookup["mapcorona"]])
        tr_entries = _extract_entries_from_map_container(idl[key_lookup["maptr"]])
        for e in cor_entries:
            e["channel"] = _extract_channel_from_id(e["id"], str(e["index"]))
        for e in tr_entries:
            e["channel"] = _extract_channel_from_id(e["id"], str(e["index"]))
        tr_cube, channels = _stack_entries_by_channel(tr_entries)
        cor_cube, channels2 = _stack_entries_by_channel(cor_entries, channel_order=channels)
        return _ViewerData(
            kind="euv",
            left_cube=tr_cube,
            right_cube=cor_cube,
            axis_kind="channel",
            axis_values=np.arange(len(channels), dtype=np.float64),
            axis_labels=channels2,
            left_label="GX (TR)",
            right_label="GX (Corona)",
            left_cmap="magma",
            right_cmap="magma",
            bunit="DN s^-1 pix^-1",
            date_obs="",
            index_header="",
        )

    raise ValueError(f"Unsupported IDL save format in {path}. Expected top-level map or mapcorona/maptr.")


def _read_render_file(path: Path) -> _ViewerData:
    ext = path.suffix.lower()
    if ext in {".h5", ".hdf5"}:
        return _read_render_h5(path)
    if ext in {".sav", ".xdr"}:
        return _read_render_sav(path)
    raise ValueError(f"Unsupported file extension: {path.suffix}")


def _build_2d_wcs_meta(index_header: str, nx: int, ny: int, date_obs: str, bunit: str) -> dict:
    meta = {
        "naxis": 2,
        "naxis1": int(nx),
        "naxis2": int(ny),
        "ctype1": "HPLN-TAN",
        "ctype2": "HPLT-TAN",
        "cunit1": "arcsec",
        "cunit2": "arcsec",
        "cdelt1": 1.0,
        "cdelt2": 1.0,
        "crpix1": (nx + 1.0) / 2.0,
        "crpix2": (ny + 1.0) / 2.0,
        "crval1": 0.0,
        "crval2": 0.0,
        "date-obs": date_obs or "",
        "bunit": bunit,
    }
    if not index_header:
        return meta
    try:
        hdr = fits.Header.fromstring(index_header, sep="\n")
    except Exception:
        return meta
    for k in ("CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CDELT1", "CDELT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "BUNIT"):
        if k in hdr:
            meta[k.lower()] = hdr[k]
    if "DATE-OBS" in hdr and not meta["date-obs"]:
        meta["date-obs"] = hdr["DATE-OBS"]
    return meta


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


def _robust_minmax(a: np.ndarray, lo_pct: float = 2.0, hi_pct: float = 98.0) -> tuple[float, float]:
    finite = np.isfinite(a)
    if not finite.any():
        return 0.0, 1.0
    vals = np.asarray(a[finite], dtype=np.float64)
    lo = float(np.nanpercentile(vals, lo_pct))
    hi = float(np.nanpercentile(vals, hi_pct))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return _safe_minmax(a)
    return lo, hi


def _nonzero_minmax(a: np.ndarray) -> tuple[float, float]:
    finite = np.isfinite(a)
    if not finite.any():
        return 0.0, 1.0
    vals = np.asarray(a[finite], dtype=np.float64)
    nz = vals[vals != 0]
    if nz.size == 0:
        return _safe_minmax(vals)
    mn = float(np.nanmin(nz))
    mx = float(np.nanmax(nz))
    if not np.isfinite(mn) or not np.isfinite(mx) or mx <= mn:
        return _safe_minmax(nz)
    return mn, mx


def _axis_display(data: _ViewerData, idx: int) -> str:
    if data["axis_kind"] == "freq":
        return f"{float(data['axis_values'][idx]):.2f} GHz"
    return f"AIA {data['axis_labels'][idx]}"


def _title(path: Path, panel_label: str, axis_text: str, date_obs: str) -> str:
    dt = f"\n{date_obs}" if date_obs else ""
    return f"{path.name} | {panel_label} @ {axis_text}{dt}"


def run_viewer(path: Path, start_index: int = 0) -> None:
    data = _read_render_file(path)
    left_cube = np.asarray(data["left_cube"], dtype=np.float64)
    right_cube = np.asarray(data["right_cube"], dtype=np.float64)
    # Viewer robustness: replace NaN/Inf pixels before plotting.
    left_cube = np.where(np.isfinite(left_cube), left_cube, 0.0)
    right_cube = np.where(np.isfinite(right_cube), right_cube, 0.0)
    if left_cube.shape != right_cube.shape:
        raise ValueError(f"Left/right cube shape mismatch: {left_cube.shape} vs {right_cube.shape}")
    if left_cube.ndim != 3:
        raise ValueError(f"Expected 3D cubes [ny, nx, n], got {left_cube.shape}")

    ny, nx, n = left_cube.shape
    base_meta = _build_2d_wcs_meta(
        index_header=str(data.get("index_header", "")),
        nx=nx,
        ny=ny,
        date_obs=str(data.get("date_obs", "")),
        bunit=str(data.get("bunit", "")),
    )
    idx0 = int(np.clip(start_index, 0, n - 1))

    left_min, left_max = _safe_minmax(left_cube)
    right_min, right_max = _safe_minmax(right_cube)

    fig = plt.figure(figsize=(12, 7.8))
    # Place log checkboxes below the range/color controls to avoid crowding the top UI.
    cb_left_ax = fig.add_axes([0.06, 0.075, 0.18, 0.05], frameon=False)
    cb_right_ax = fig.add_axes([0.58, 0.075, 0.18, 0.05], frameon=False)
    cb_left_ax.set_xticks([])
    cb_left_ax.set_yticks([])
    cb_right_ax.set_xticks([])
    cb_right_ax.set_yticks([])
    ax_slider = fig.add_axes([0.24, 0.94, 0.62, 0.03])

    left_data = left_cube[:, :, idx0]
    right_data = right_cube[:, :, idx0]

    left_meta = dict(base_meta)
    right_meta = dict(base_meta)
    left_meta["content"] = str(data["left_label"])
    right_meta["content"] = str(data["right_label"])
    if data["axis_kind"] == "freq":
        left_meta["freqghz"] = float(data["axis_values"][idx0])
        right_meta["freqghz"] = float(data["axis_values"][idx0])
    else:
        left_meta["channel"] = str(data["axis_labels"][idx0])
        right_meta["channel"] = str(data["axis_labels"][idx0])

    m_left0 = sunpy.map.Map(left_data, left_meta)
    m_right0 = sunpy.map.Map(right_data, right_meta)
    ax_left = fig.add_axes([0.06, 0.27, 0.38, 0.60], projection=m_left0)
    ax_right = fig.add_axes([0.56, 0.27, 0.38, 0.60], projection=m_right0)
    ax_left_rng = fig.add_axes([0.06, 0.14, 0.36, 0.035])
    ax_right_rng = fig.add_axes([0.58, 0.14, 0.36, 0.035])
    left_info = fig.text(0.06, 0.182, "", ha="left", va="bottom")
    right_info = fig.text(0.58, 0.182, "", ha="left", va="bottom")

    left_init_min, left_init_max = _robust_minmax(left_data)
    right_init_min, right_init_max = _robust_minmax(right_data)
    left_nz_min, left_nz_max = _nonzero_minmax(left_data)
    right_nz_min, right_nz_max = _nonzero_minmax(right_data)

    im_left = m_left0.plot(
        axes=ax_left,
        cmap=str(data["left_cmap"]),
        norm=_norm_for_data(left_data, left_init_min, left_init_max, log=False),
        interpolation="nearest",
    )
    im_right = m_right0.plot(
        axes=ax_right,
        cmap=str(data["right_cmap"]),
        norm=_norm_for_data(right_data, right_init_min, right_init_max, log=False),
        interpolation="nearest",
    )

    axis_text = _axis_display(data, idx0)
    ax_left.set_title(_title(path, str(data["left_label"]), axis_text, str(data.get("date_obs", ""))))
    ax_right.set_title(_title(path, str(data["right_label"]), axis_text, str(data.get("date_obs", ""))))
    ax_left.set_xlabel("Solar X [arcsec]")
    ax_left.set_ylabel("Solar Y [arcsec]")
    ax_right.set_xlabel("Solar X [arcsec]")
    ax_right.set_ylabel("Solar Y [arcsec]")

    cbar_left = fig.colorbar(im_left, ax=ax_left)
    cbar_right = fig.colorbar(im_right, ax=ax_right)
    cbar_left.set_label(str(data["bunit"]))
    cbar_right.set_label(str(data["bunit"]))

    slider_label = "Frequency index" if data["axis_kind"] == "freq" else "Channel index"
    axis_slider = Slider(ax=ax_slider, label=slider_label, valmin=0, valmax=n - 1, valinit=idx0, valstep=1)
    left_range = RangeSlider(
        ax=ax_left_rng, label="", valmin=left_nz_min, valmax=left_nz_max, valinit=(left_nz_min, left_nz_max)
    )
    right_range = RangeSlider(
        ax=ax_right_rng, label="", valmin=right_nz_min, valmax=right_nz_max, valinit=(right_nz_min, right_nz_max)
    )
    left_range.valtext.set_visible(False)
    right_range.valtext.set_visible(False)

    cb_left = CheckButtons(cb_left_ax, [f"{data['left_label']} log"], [False])
    cb_right = CheckButtons(cb_right_ax, [f"{data['right_label']} log"], [False])
    state = {"idx": idx0, "left_log": False, "right_log": False, "updating_ranges": False}

    def _set_slider_bounds(slider: RangeSlider, vmin: float, vmax: float) -> None:
        # Keep a valid span for constant maps.
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            vmin, vmax = 0.0, 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        slider.valmin = float(vmin)
        slider.valmax = float(vmax)
        slider.ax.set_xlim(float(vmin), float(vmax))

    def _reset_ranges_for_idx(idx: int) -> None:
        state["updating_ranges"] = True
        try:
            lmn, lmx = _nonzero_minmax(left_cube[:, :, idx])
            rmn, rmx = _nonzero_minmax(right_cube[:, :, idx])
            _set_slider_bounds(left_range, lmn, lmx)
            _set_slider_bounds(right_range, rmn, rmx)
            left_range.set_val((lmn, lmx))
            right_range.set_val((rmn, rmx))
        finally:
            state["updating_ranges"] = False

    def redraw() -> None:
        idx = state["idx"]
        left = left_cube[:, :, idx]
        right = right_cube[:, :, idx]
        im_left.set_data(left)
        im_right.set_data(right)

        left_vmin, left_vmax = left_range.val
        right_vmin, right_vmax = right_range.val
        im_left.set_norm(_norm_for_data(left, float(left_vmin), float(left_vmax), bool(state["left_log"])))
        im_right.set_norm(_norm_for_data(right, float(right_vmin), float(right_vmax), bool(state["right_log"])))
        left_info.set_text(f"{data['left_label']} range [{data['bunit']}]: ({left_vmin:.3g}, {left_vmax:.3g})")
        right_info.set_text(f"{data['right_label']} range [{data['bunit']}]: ({right_vmin:.3g}, {right_vmax:.3g})")

        axis_text = _axis_display(data, idx)
        ax_left.set_title(_title(path, str(data["left_label"]), axis_text, str(data.get("date_obs", ""))))
        ax_right.set_title(_title(path, str(data["right_label"]), axis_text, str(data.get("date_obs", ""))))
        fig.canvas.draw_idle()

    def _on_axis(v):
        state["idx"] = int(v)
        _reset_ranges_for_idx(state["idx"])
        redraw()

    def _on_left_range(_):
        if not state["updating_ranges"]:
            redraw()

    def _on_right_range(_):
        if not state["updating_ranges"]:
            redraw()

    axis_slider.on_changed(_on_axis)
    left_range.on_changed(_on_left_range)
    right_range.on_changed(_on_right_range)
    cb_left.on_clicked(lambda _: (state.__setitem__("left_log", not state["left_log"]), redraw()))
    cb_right.on_clicked(lambda _: (state.__setitem__("right_log", not state["right_log"]), redraw()))
    redraw()
    plt.show()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Interactive viewer for rendered MW/EUV map files (.h5 or IDL .sav).")
    p.add_argument("map_path", type=Path, help="Path to rendered map file (.h5/.hdf5/.sav/.xdr).")
    p.add_argument("--start-index", type=int, default=0, help="Initial frequency/channel index.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_viewer(args.map_path.expanduser().resolve(), start_index=args.start_index)


if __name__ == "__main__":
    main()
