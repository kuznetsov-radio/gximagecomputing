from __future__ import annotations

import numpy as np


def gx_voxelid(
    test_id: int | np.ndarray | None = None,
    *,
    chromo: bool = False,
    tr: bool = False,
    corona: bool = False,
    euvtr: bool = False,
    in_: bool = False,
    nw: bool = False,
    enw: bool = False,
    fa: bool = False,
    pl: bool = False,
    pen: bool = False,
    umb: bool = False,
    tube: bool = False,
    layer: bool = False,
    mask: bool = False,
) -> np.uint32 | np.ndarray:
    """IDL-equivalent voxel bitmask builder/tester.

    Port of ``gx_voxelid.pro`` behavior.
    """
    vid = np.uint32(0)

    if chromo:
        vid += np.uint32(1)
    if tr:
        vid += np.uint32(2)
    if corona:
        vid += np.uint32(4)
    if euvtr:
        vid += np.uint32(8)

    if tube:
        vid += np.uint32(255 << 8)
    if layer:
        vid += np.uint32(255 << 16)
    if mask:
        vid += np.uint32(255 << 24)

    if in_:
        vid += np.uint32(1 << 24)
    if nw:
        vid += np.uint32(2 << 24)
    if enw:
        vid += np.uint32(3 << 24)
    if fa:
        vid += np.uint32(4 << 24)
    if pl:
        vid += np.uint32(5 << 24)
    if pen:
        vid += np.uint32(6 << 24)
    if umb:
        vid += np.uint32(7 << 24)

    if test_id is None:
        return vid
    return vid & np.asarray(test_id, dtype=np.uint32)


def _has_tag(box, name: str) -> bool:
    if isinstance(box, dict):
        return name in box
    return hasattr(box, name)


def _get_tag(box, name: str, default=None):
    if isinstance(box, dict):
        return box.get(name, default)
    return getattr(box, name, default)


def _find_field(box, candidates: tuple[str, ...]):
    for name in candidates:
        if _has_tag(box, name):
            return _get_tag(box, name)
    return None


def _gx_tr_height_sunradius() -> float:
    # IDL parity for gx_tr_height(): 2000 / gx_rsun(km)
    return 2000.0 / 696000.0


def _normalize_box_for_id(box):
    """Accept either flat dict or grouped dict with 'corona'/'chromo'."""
    if not isinstance(box, dict):
        return box
    if "corona" not in box and "chromo" not in box:
        return box
    chromo = box.get("chromo", {}) if isinstance(box.get("chromo", {}), dict) else {}
    corona = box.get("corona", {}) if isinstance(box.get("corona", {}), dict) else {}
    merged = {}
    for key in ("bx", "by", "bz", "Bx", "By", "Bz", "bcube", "dr", "start_idx", "startidx"):
        if key in corona:
            merged[key] = corona[key]
    for key in ("chromo_idx", "chromo_t", "chromo_n", "chromo_layers", "corona_base"):
        if key in chromo:
            merged[key] = chromo[key]
    for key in (
        "bx",
        "by",
        "bz",
        "Bx",
        "By",
        "Bz",
        "bcube",
        "dr",
        "start_idx",
        "startidx",
        "chromo_idx",
        "chromo_t",
        "chromo_n",
        "chromo_layers",
        "corona_base",
    ):
        if key in box and key not in merged:
            merged[key] = box[key]
    return merged if merged else box


def gx_box2id(box, tr_mask: np.ndarray | None = None):
    """Build voxel_id 3D array from a box structure (IDL parity)."""
    if box is None:
        return None
    box = _normalize_box_for_id(box)

    bx = _find_field(box, ("Bx", "bx"))
    by = _find_field(box, ("By", "by"))
    bz = _find_field(box, ("Bz", "bz"))
    valid_bxbybz = bx is not None and by is not None and bz is not None
    bcube = _find_field(box, ("bcube",))
    valid_b = bcube is not None
    if not (valid_bxbybz or valid_b):
        return None

    dr = np.asarray(_find_field(box, ("dr",)), dtype=float)
    if dr.size < 3:
        raise ValueError("box.dr must contain at least 3 elements.")

    if _has_tag(box, "corona_base"):
        corona_base = int(_get_tag(box, "corona_base"))
    else:
        start_idx = _find_field(box, ("startidx", "start_idx"))
        if start_idx is not None:
            start_idx = np.asarray(start_idx)
            nonzero_vals = start_idx[start_idx != 0]
            if nonzero_vals.size > 0:
                _, _, z_idx = np.unravel_index(nonzero_vals.astype(np.int64), start_idx.shape, order="F")
                corona_base = int(np.min(z_idx))
            else:
                corona_base = int(np.ceil(_gx_tr_height_sunradius() / dr[2]))
        else:
            corona_base = int(np.ceil(_gx_tr_height_sunradius() / dr[2]))

    chromo_layers = int(_get_tag(box, "chromo_layers", corona_base))

    if valid_bxbybz:
        sz = np.asarray(np.shape(bx), dtype=int)
    else:
        sz = np.asarray(np.shape(bcube[..., 0]), dtype=int)
    if sz.size != 3:
        raise ValueError("Expected magnetic field volume with shape (nx, ny, nz).")

    out_sz = sz.copy()
    if corona_base > 0:
        out_sz[2] = out_sz[2] - corona_base + chromo_layers
    nx, ny, nz = (int(out_sz[0]), int(out_sz[1]), int(out_sz[2]))

    vid = np.zeros((nx, ny, nz), dtype=np.uint32)
    tr = np.zeros((nx, ny), dtype=np.uint32)

    has_combo = _has_tag(box, "chromo_idx") and _has_tag(box, "chromo_t") and _has_tag(box, "chromo_n")
    if has_combo:
        vol_t = np.zeros((nx, ny, nz), dtype=float)
        vol_n = np.zeros((nx, ny, nz), dtype=float)

        chromo_idx = np.asarray(_get_tag(box, "chromo_idx"), dtype=np.int64)
        chromo_t = np.asarray(_get_tag(box, "chromo_t"), dtype=float)
        chromo_n = np.asarray(_get_tag(box, "chromo_n"), dtype=float)

        vtf = vol_t.ravel(order="F")
        vnf = vol_n.ravel(order="F")
        nfill = min(chromo_idx.size, chromo_t.size, chromo_n.size)
        valid = (chromo_idx[:nfill] >= 0) & (chromo_idx[:nfill] < vtf.size)
        idx = chromo_idx[:nfill][valid]
        vtf[idx] = chromo_t[:nfill][valid]
        vnf[idx] = chromo_n[:nfill][valid]
        vol_t = vtf.reshape(vol_t.shape, order="F")
        vol_n = vnf.reshape(vol_n.shape, order="F")

        for i in range(nx):
            for j in range(ny):
                k = np.where((vol_n[i, j, :] != 0) & (vol_t[i, j, :] != 0))[0]
                tr[i, j] = np.uint32((int(np.max(k)) + 1) if k.size else 0)
    else:
        tr[:, :] = np.uint32(max(corona_base, 0))

    max_tr = int(np.max(tr))
    if max_tr >= 0:
        k_end = min(max_tr + 1, nz)
        vid[:, :, 0:k_end] = gx_voxelid(chromo=True)

    tr_euv = gx_voxelid(tr=True, euvtr=True)
    cor = gx_voxelid(corona=True)
    for i in range(nx):
        for j in range(ny):
            k = int(tr[i, j])
            if 0 <= k < nz:
                vid[i, j, k] = np.uint32(vid[i, j, k] + tr_euv)
                vid[i, j, k:nz] = np.uint32(vid[i, j, k:nz] + cor)

    cor_chromo_id = gx_voxelid(chromo=True, corona=True)
    coridx = np.where(vid.ravel(order="F") == cor_chromo_id)[0]
    if coridx.size > 0:
        flat = vid.ravel(order="F")
        valid = (coridx > 0) & (coridx < flat.size - 1)
        coridx_v = coridx[valid]
        if coridx_v.size > 0:
            edge = np.where((flat[coridx_v + 1] == 1) | (flat[coridx_v - 1] == 1))[0]
            if edge.size > 0:
                flat[coridx_v[edge]] = np.uint32(flat[coridx_v[edge]] + tr_euv)
                vid = flat.reshape(vid.shape, order="F")

    if tr_mask is not None:
        tr_mask = np.asarray(tr_mask)
        if tr_mask.shape[0] == nx and tr_mask.shape[1] == ny:
            mask = np.uint32(1) - np.asarray(tr_mask, dtype=np.uint32)
            tr_euv_mask = gx_voxelid(tr=True, euvtr=True)
            flat = vid.ravel(order="F")
            tr_idx = np.where((flat & tr_euv_mask) != 0)[0]
            if tr_idx.size > 0:
                xi, yi, _ = np.unravel_index(tr_idx, vid.shape, order="F")
                euv = gx_voxelid(euvtr=True)
                flat[tr_idx] = np.uint32(flat[tr_idx] - mask[xi, yi] * euv)
                vid = flat.reshape(vid.shape, order="F")

    return vid
