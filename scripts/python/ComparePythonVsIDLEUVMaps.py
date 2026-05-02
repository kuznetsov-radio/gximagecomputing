#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import ComparePythonVsIDLMaps as base


DEFAULT_OUTDIR = Path(tempfile.gettempdir()) / "gximagecomputing_validation_groundtruth"

# Repoint defaults to EUV artifacts while reusing the shared comparison implementation.
base.DEFAULT_PY_H5 = DEFAULT_OUTDIR / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5_py_euv_maps.h5"
base.DEFAULT_IDL_SAV = DEFAULT_OUTDIR / "hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.sav_idl_euv_maps.sav"


def main() -> None:
    argv = sys.argv[1:]
    if "--kind" not in argv:
        sys.argv = [sys.argv[0], "--kind", "euv", *argv]
    base.main()


if __name__ == "__main__":
    main()
