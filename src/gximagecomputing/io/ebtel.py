from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import scipy.io as sio


# Keep EBTEL resolution explicit and portable:
# - CLI/API explicit path takes precedence.
# - Optional environment override:
#     GXIMAGECOMPUTING_EBTEL_PATH=/path/to/ebtel_ss.sav
EBTEL_CANDIDATES: list[Path] = []


def resolve_ebtel_path(explicit_path: Path | None = None) -> Path:
    candidates = [explicit_path] if explicit_path is not None else []
    env_path = os.environ.get("GXIMAGECOMPUTING_EBTEL_PATH")
    if env_path:
        candidates.append(Path(env_path))
    candidates.extend(EBTEL_CANDIDATES)
    tried = []
    for path in candidates:
        if path is None:
            continue
        p = Path(path).expanduser()
        tried.append(str(p))
        if p.exists():
            return p
    tried_msg = "\n".join(f"- {p}" for p in tried)
    raise FileNotFoundError(f"No EBTEL file found. Paths checked:\n{tried_msg}")


def decode_if_bytes(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8")
    return value


def load_ebtel(ebtel_file: str | os.PathLike[str]):
    ebtel_data = sio.readsav(str(ebtel_file))
    s = ebtel_data["lrun"].shape

    ebtel_dtypes = [
        ("DEM_on", np.int32),
        ("DDM_on", np.int32),
        ("NQ", np.int32),
        ("NL", np.int32),
        ("NT", np.int32),
        ("Qrun", np.float32, ebtel_data["qrun"].shape),
        ("Lrun", np.float32, ebtel_data["lrun"].shape),
        ("logtdem", np.float32, ebtel_data["logtdem"].size),
    ]
    if "dem_cor_run" in ebtel_data.keys():
        ebtel_dtypes.append(("DEM_cor_run", np.float32, ebtel_data["dem_cor_run"].shape))
    if "ddm_cor_run" in ebtel_data.keys():
        ebtel_dtypes.append(("DDM_cor_run", np.float32, ebtel_data["ddm_cor_run"].shape))
    if "dem_tr_run" in ebtel_data.keys():
        ebtel_dtypes.append(("DEM_TR_RUN", np.float32, ebtel_data["dem_tr_run"].shape))
    if "ddm_tr_run" in ebtel_data.keys():
        ebtel_dtypes.append(("DDM_TR_RUN", np.float32, ebtel_data["ddm_tr_run"].shape))

    ebtel_dt = np.dtype(ebtel_dtypes)
    ebtel_c = np.zeros(1, dtype=ebtel_dt)

    ebtel_c["DEM_on"] = 1 if "dem_cor_run" in ebtel_data.keys() else 0
    ebtel_c["DDM_on"] = 1 if "ddm_cor_run" in ebtel_data.keys() else 0
    # Python/C ABI in this project expects NQ/NL in this order.
    ebtel_c["NQ"] = s[1]
    ebtel_c["NL"] = s[0]
    ebtel_c["NT"] = ebtel_data["logtdem"].size
    ebtel_c["Qrun"] = np.float32(ebtel_data["qrun"])
    ebtel_c["Lrun"] = np.float32(ebtel_data["lrun"])
    ebtel_c["logtdem"] = np.float32(ebtel_data["logtdem"])

    if "dem_cor_run" in ebtel_data.keys():
        ebtel_c["DEM_cor_run"] = np.float32(ebtel_data["DEM_cor_run"])
    if "ddm_cor_run" in ebtel_data.keys():
        ebtel_c["DDM_cor_run"] = np.float32(ebtel_data["DDM_cor_run"])
    if "dem_tr_run" in ebtel_data.keys():
        ebtel_c["DEM_TR_RUN"] = np.float32(ebtel_data["DEM_TR_RUN"])
    if "ddm_tr_run" in ebtel_data.keys():
        ebtel_c["DDM_TR_RUN"] = np.float32(ebtel_data["DDM_TR_RUN"])

    return ebtel_c, ebtel_dt


def load_ebtel_none():
    """Return a minimal EBTEL struct with DEM/DDM disabled."""
    ebtel_dt = np.dtype(
        [
            ("DEM_on", np.int32),
            ("DDM_on", np.int32),
            ("NQ", np.int32),
            ("NL", np.int32),
            ("NT", np.int32),
            ("Qrun", np.float32, (1, 1)),
            ("Lrun", np.float32, (1, 1)),
            ("logtdem", np.float32, (1,)),
        ]
    )
    ebtel_c = np.zeros(1, dtype=ebtel_dt)
    ebtel_c["DEM_on"] = 0
    ebtel_c["DDM_on"] = 0
    ebtel_c["NQ"] = 1
    ebtel_c["NL"] = 1
    ebtel_c["NT"] = 1
    return ebtel_c, ebtel_dt
