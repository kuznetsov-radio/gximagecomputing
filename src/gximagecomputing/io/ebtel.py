from __future__ import annotations

import os
from pathlib import Path

import numpy as np


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
