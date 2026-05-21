from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
from astropy.io import fits

from gxrender.workflows.render_euv import _preview_euv


def _decode_scalar(value) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _read_meta_text(meta: h5py.Group | None, key: str) -> str:
    if meta is None or key not in meta:
        return ""
    return _decode_scalar(meta[key][()])


def _build_header(meta: h5py.Group | None, nx: int, ny: int) -> fits.Header:
    wcs_header = _read_meta_text(meta, "wcs_header")
    index_header = _read_meta_text(meta, "index_header")
    header_text = wcs_header or index_header
    if header_text:
        try:
            return fits.Header.fromstring(header_text, sep="\n")
        except Exception:
            pass

    hdr = fits.Header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = int(nx)
    hdr["NAXIS2"] = int(ny)
    hdr["CTYPE1"] = "HPLN-TAN"
    hdr["CTYPE2"] = "HPLT-TAN"
    hdr["CUNIT1"] = "arcsec"
    hdr["CUNIT2"] = "arcsec"
    hdr["CRPIX1"] = (nx + 1.0) / 2.0
    hdr["CRPIX2"] = (ny + 1.0) / 2.0
    if meta is not None:
        def _meta_float(name: str, default: float) -> float:
            if name in meta:
                try:
                    return float(meta[name][()])
                except Exception:
                    return default
            return default

        hdr["CRVAL1"] = _meta_float("xc_arcsec", 0.0)
        hdr["CRVAL2"] = _meta_float("yc_arcsec", 0.0)
        hdr["CDELT1"] = _meta_float("dx_arcsec", 1.0)
        hdr["CDELT2"] = _meta_float("dy_arcsec", 1.0)
        date_obs = _read_meta_text(meta, "date_obs")
        if date_obs:
            hdr["DATE-OBS"] = date_obs
        bunit = _read_meta_text(meta, "bunit") or "DN s^-1 pix^-1"
        hdr["BUNIT"] = bunit
    return hdr


def _resolve_component_indices(component_ids: list[str]) -> tuple[int, int]:
    upper = [value.upper() for value in component_ids]
    cor_idx = next((i for i, value in enumerate(upper) if value.startswith("COR")), 0)
    tr_idx = next((i for i, value in enumerate(upper) if value.startswith("TR")), 1 if len(upper) > 1 else 0)
    return cor_idx, tr_idx


def _resolve_channel_index(channels: list[str], channel_id: str | None, channel_index: int) -> int:
    if not channels:
        raise ValueError("No channels found in rendered map file.")
    if channel_id:
        wanted = str(channel_id).strip().upper()
        if wanted in {channel.upper() for channel in channels}:
            for idx, channel in enumerate(channels):
                if channel.upper() == wanted:
                    return idx
        if wanted.startswith("A") and wanted[1:] in {channel.lstrip("A").upper() for channel in channels}:
            raw = wanted[1:]
            for idx, channel in enumerate(channels):
                if channel.lstrip("A").upper() == raw:
                    return idx
        raise ValueError(f"Requested channel {channel_id!r} not found. Available channels: {','.join(channels)}")
    return int(np.clip(int(channel_index), 0, len(channels) - 1))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Regenerate an EUV preview PNG from a rendered H5 map file.")
    parser.add_argument("map_path", type=Path, help="Rendered EUV map H5 path.")
    parser.add_argument("--channel-index", type=int, default=0, help="Channel index to preview (default: 0).")
    parser.add_argument("--channel-id", type=str, default=None, help="Explicit channel label (for example A171).")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output PNG path. Default: <map-stem>_preview_<channel>.png beside the map.",
    )
    parser.add_argument("--show", action="store_true", help="Display preview in-memory with matplotlib.")
    parser.add_argument("--no-save", action="store_true", help="Do not save PNG to disk.")
    parser.add_argument("--title", type=str, default=None, help="Optional preview title override.")
    parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Apply log10 scaling to preview PNG for enhanced dynamic range visualization.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    map_path = args.map_path.expanduser().resolve()
    if not map_path.exists():
        raise FileNotFoundError(f"Map file not found: {map_path}")
    if args.no_save and args.output is not None:
        raise ValueError("--output cannot be used together with --no-save")

    with h5py.File(map_path, "r") as handle:
        if "maps" not in handle or "data" not in handle["maps"]:
            raise ValueError(f"{map_path} is not a rendered map file with maps/data")
        maps = handle["maps"]
        data = np.asarray(maps["data"], dtype=np.float64)
        if data.ndim != 4 or data.shape[-1] != 2:
            raise ValueError(f"Expected data shape [nx, ny, n, 2], got {data.shape}")
        if "channel_ids" not in maps:
            raise ValueError("This file does not contain EUV channel_ids; expected EUV map H5 artifact.")

        channels = [_decode_scalar(value) for value in np.asarray(maps["channel_ids"]) ]
        component_ids = [_decode_scalar(value) for value in np.asarray(maps.get("component_ids", np.array(["CORONA", "TR"])))]
        channel_idx = _resolve_channel_index(channels, args.channel_id, args.channel_index)

        disp = np.transpose(data, (1, 0, 2, 3))
        cor_idx, tr_idx = _resolve_component_indices(component_ids)
        flux_cor = disp[:, :, channel_idx, cor_idx]
        flux_tr = disp[:, :, channel_idx, tr_idx]

        meta = handle.get("metadata")
        ny, nx = flux_cor.shape
        header = _build_header(meta, nx=nx, ny=ny)
        instrument = _read_meta_text(meta, "instrument") or "EUV"
        observer_name = _read_meta_text(meta, "observer_name")

    channel = channels[channel_idx]
    preview_path = None
    if not args.no_save:
        preview_path = args.output.expanduser().resolve() if args.output is not None else map_path.with_name(
            f"{map_path.stem}_preview_{channel}.png"
        )

    _preview_euv(
        flux_cor=flux_cor,
        flux_tr=flux_tr,
        instrument=instrument,
        channel=channel,
        out_png=preview_path,
        title=str(args.title or map_path.stem),
        base_wcs_header=header,
        observer_text=observer_name,
        show=bool(args.show),
        log_scale=bool(args.log_scale),
    )

    if preview_path is not None:
        print(f"preview_png: {preview_path}")
    else:
        print(f"preview_in_memory: instrument={instrument} channel={channel}")


if __name__ == "__main__":
    main()
