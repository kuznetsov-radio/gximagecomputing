#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import tempfile
from pathlib import Path

import h5py
import numpy as np
from scipy.io import readsav
from gximagecomputing.io.voxel_id import gx_box2id


def _decode_if_bytes(v):
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8")
    return v


def _as_scalar(v):
    arr = np.asarray(v)
    if arr.shape == ():
        return arr.item()
    return arr.flat[0]


def _load_box(sav_path: Path):
    data = readsav(str(sav_path), verbose=False)
    if "box" not in data:
        raise ValueError("Input SAV does not contain 'box'.")
    return data["box"].flat[0]


def _ensure_group(f: h5py.File, name: str):
    if name in f:
        return f[name]
    return f.create_group(name)


def build_h5_from_sav(template_h5: Path, sav_path: Path, out_h5: Path) -> None:
    shutil.copy2(template_h5, out_h5)
    box = _load_box(sav_path)
    index = box["INDEX"][0] if "INDEX" in box.dtype.names else None

    # Ground-truth mapped arrays from SAV (IDL conventions).
    dr = np.asarray(box["DR"], dtype=np.float64)
    dz = np.asarray(box["DZ"], dtype=np.float64).transpose((1, 2, 0))
    bcube = np.asarray(box["BCUBE"], dtype=np.float32).transpose((2, 3, 1, 0))
    chromo_bcube = np.asarray(box["CHROMO_BCUBE"], dtype=np.float32).transpose((2, 3, 1, 0))

    av_field = np.asarray(box["AVFIELD"], dtype=np.float64).T.reshape(-1, order="F")
    phys_length = np.asarray(box["PHYSLENGTH"], dtype=np.float64).T.reshape(-1, order="F")
    voxel_status = np.asarray(box["STATUS"], dtype=np.uint8).T.reshape(-1, order="F")
    start_idx = np.asarray(box["STARTIDX"], dtype=np.int32).T.reshape(-1, order="F")
    end_idx = np.asarray(box["ENDIDX"], dtype=np.int32).T.reshape(-1, order="F")

    chromo_idx = np.asarray(box["CHROMO_IDX"], dtype=np.int64)
    chromo_n = np.asarray(box["CHROMO_N"], dtype=np.float32)
    chromo_t = np.asarray(box["CHROMO_T"], dtype=np.float32)
    n_p = np.asarray(box["N_P"], dtype=np.float32)
    n_hi = np.asarray(box["N_HI"], dtype=np.float32)
    n_htot = np.asarray(box["N_HTOT"], dtype=np.float32)
    tr = np.asarray(box["TR"], dtype=np.int64)
    tr_h = np.asarray(box["TR_H"], dtype=np.float32)
    chromo_layers = int(box["CHROMO_LAYERS"])
    corona_base = int(box["CORONA_BASE"])

    base = box["BASE"][0]
    base_bx = np.asarray(base["BX"], dtype=np.float64)
    base_by = np.asarray(base["BY"], dtype=np.float64)
    base_bz = np.asarray(base["BZ"], dtype=np.float64)
    base_ic = np.asarray(base["IC"], dtype=np.float64)
    base_chromo_mask = np.asarray(base["CHROMO_MASK"], dtype=np.int32)

    with h5py.File(out_h5, "r+") as f:
        g_chromo = _ensure_group(f, "chromo")
        g_base = _ensure_group(f, "base")
        g_grid = _ensure_group(f, "grid")
        g_meta = _ensure_group(f, "metadata")

        # Write/overwrite key chromo datasets.
        write_map = {
            "dr": dr,
            "dz": dz,
            "bcube": bcube,
            "chromo_bcube": chromo_bcube,
            "av_field": av_field,
            "phys_length": phys_length,
            "voxel_status": voxel_status,
            "start_idx": start_idx,
            "end_idx": end_idx,
            "chromo_idx": chromo_idx,
            "chromo_n": chromo_n,
            "chromo_t": chromo_t,
            "n_p": n_p,
            "n_hi": n_hi,
            "n_htot": n_htot,
            "tr": tr,
            "tr_h": tr_h,
            "chromo_layers": np.array(chromo_layers, dtype=np.int64),
            "corona_base": np.array(corona_base, dtype=np.int64),
            "chromo_mask": base_chromo_mask,
        }
        for key, arr in write_map.items():
            if key in g_chromo:
                del g_chromo[key]
            g_chromo.create_dataset(key, data=arr)

        # Ensure HDF loader uses SAV source-of-truth pointing/time metadata.
        if index is not None:
            g_chromo.attrs["lon"] = float(_as_scalar(index["CRVAL1"]))
            g_chromo.attrs["lat"] = float(_as_scalar(index["CRVAL2"]))
            g_chromo.attrs["dsun_obs"] = float(_as_scalar(index["DSUN_OBS"]))
            g_chromo.attrs["obs_time"] = _decode_if_bytes(_as_scalar(index["DATE_OBS"]))

        # Base maps from SAV source of truth.
        for key, arr in {
            "bx": base_bx,
            "by": base_by,
            "bz": base_bz,
            "ic": base_ic,
            "chromo_mask": base_chromo_mask,
        }.items():
            if key in g_base:
                del g_base[key]
            g_base.create_dataset(key, data=arr)

        # Build voxel_id on the fly from source-of-truth arrays.
        flat_box = {
            "bcube": bcube,
            "dr": dr,
            "start_idx": start_idx,
            "chromo_idx": chromo_idx,
            "chromo_t": chromo_t,
            "chromo_n": chromo_n,
            "chromo_layers": chromo_layers,
            "corona_base": corona_base,
        }
        voxel_id = gx_box2id(flat_box)
        if voxel_id is None:
            raise RuntimeError("gx_box2id returned None; could not compute voxel_id.")

        if "voxel_id" in g_grid:
            del g_grid["voxel_id"]
        g_grid.create_dataset("voxel_id", data=np.asarray(voxel_id, dtype=np.uint32))

        if "dx" in g_grid:
            del g_grid["dx"]
        if "dy" in g_grid:
            del g_grid["dy"]
        if "dz" in g_grid:
            del g_grid["dz"]
        g_grid.create_dataset("dx", data=np.float64(dr[0]))
        g_grid.create_dataset("dy", data=np.float64(dr[1]))
        g_grid.create_dataset("dz", data=dz)

        # Metadata update for traceability.
        execute = _decode_if_bytes(box["EXECUTE"]) if "EXECUTE" in box.dtype.names else ""
        model_id = _decode_if_bytes(box["ID"]) if "ID" in box.dtype.names else out_h5.stem
        for key, val in {
            "execute": f"{execute} [sav->h5 ground truth copy]",
            "id": model_id,
            "disambiguation": "IDL",
        }.items():
            if key in g_meta:
                del g_meta[key]
            g_meta.create_dataset(key, data=np.bytes_(str(val)))


def parse_args():
    p = argparse.ArgumentParser(description="Build pyAMPP-like H5 from CHR SAV as source-of-truth.")
    p.add_argument(
        "--template-h5",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.h5",
        help="Existing H5 template used for non-CHR groups/layout.",
    )
    p.add_argument(
        "--sav-path",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "test_data" / "test.chr.sav",
    )
    p.add_argument(
        "--out-h5",
        type=Path,
        default=Path(tempfile.gettempdir()) / "chr_from_sav_groundtruth.h5",
    )
    return p.parse_args()


def main():
    args = parse_args()
    args.out_h5.parent.mkdir(parents=True, exist_ok=True)
    build_h5_from_sav(args.template_h5, args.sav_path, args.out_h5)
    print("Outputs:")
    print(f"- out_h5: {args.out_h5}")


if __name__ == "__main__":
    main()
