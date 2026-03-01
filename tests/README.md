# Tests / Internal Validation Workflows

This document is for repository-internal regression and parity workflows.
These procedures are useful for development/debugging and are not part of the
public user-facing quick-start documentation.

## Contents

- SAV↔H5 model parity checks
- IDL/Python MW/EUV map parity comparisons
- ComputeEUV pre-DLL input dump/compare workflows

## One-Command SAV↔H5 Parity Check

Run strict regression parity (rebuild H5 from SAV, then compare loader outputs field-by-field):

```bash
make parity-roundtrip
```

Defaults:
- `SAV_PATH=test_data/test.chr.sav`
- `H5_PATH=test_data/test.chr.h5`
- `ATOL=0`, `RTOL=0`

Override example:

```bash
make parity-roundtrip \
  SAV_PATH=/path/to/model.NAS.CHR.sav \
  H5_PATH=/tmp/model.NAS.CHR.h5
```

If your `test_data` is under Dropbox/iCloud and file locking interferes, prefer a temporary output file:

```bash
make parity-roundtrip H5_PATH=/tmp/test.chr.h5
```

## IDL/Python EUV Parity Mode (Same Input = Same Output)

For strict parity testing between IDL and Python, use identical:

- model file
- response function (`RESPonsefile` / `--response-sav`)
- EBTEL table
- map geometry (`XC/YC/DX/DY/NX/NY` or Python CLI equivalents)
- observer metadata (`DSUN/LONC/B0SUN` in IDL, `--dsun-cm/--lonc-deg/--b0sun-deg` in Python)

The IDL examples accept:

- `RESPonsefile=`
- `DSUN=`
- `LONC=`
- `B0SUN=`

and the Python MW/EUV CLIs accept the matching observer override flags.

### Example EUV parity commands (forced observer + geometry)

IDL:

```idl
RenderExampleEUV, $
  MODelfile='/path/to/model.chr.sav', $
  EBTELfile='/path/to/ebtel.sav', $
  RESPonsefile='/path/to/aia_response.sav', $
  DSUN=14763359700479.328d, $
  LONC=-17.0574058213d, $
  B0SUN=1.4406505929155138d, $
  XC=-279.97540414889585d, YC=-229.983277241489d, $
  OUTfile='/tmp/gximagecomputing_validation_groundtruth/idl_euv_maps_forced.sav', $
  /NO_PLOT
```

Python:

```bash
SUNPY_CONFIGDIR=/tmp/.sunpy-config MPLCONFIGDIR=/tmp/.mplconfig XDG_CACHE_HOME=/tmp/.cache \
PYTHONPATH=src python examples/python/cli/RenderExampleEUV.py \
  --model-path /path/to/model.chr.sav \
  --model-format sav \
  --ebtel-path /path/to/ebtel.sav \
  --response-sav /path/to/aia_response.sav \
  --dsun-cm 14763359700479.328 \
  --lonc-deg -17.0574058213 \
  --b0sun-deg 1.4406505929155138 \
  --xc -279.97540414889585 \
  --yc -229.983277241489 \
  --output-dir /tmp/gximagecomputing_validation_groundtruth
```

Compare EUV outputs (correct `CORONA/TR` labels):

```bash
SUNPY_CONFIGDIR=/tmp/.sunpy-config MPLCONFIGDIR=/tmp/.mplconfig XDG_CACHE_HOME=/tmp/.cache \
PYTHONPATH=src python tests/ComparePythonVsIDLEUVMaps.py \
  --python-h5 /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_py_euv_maps.h5 \
  --idl-sav /tmp/gximagecomputing_validation_groundtruth/idl_euv_maps_forced.sav
```

## Compare Scripts

MW/EUV map comparisons:
- `tests/ComparePythonVsIDLMaps.py` (shared implementation; supports `--kind auto|mw|euv`)
- `tests/ComparePythonVsIDLEUVMaps.py` (EUV wrapper with `--kind euv` default)

ComputeEUV pre-DLL input parity:
- `tests/DumpComputeEUVInputs.py` (Python-side input dump)
- `tests/dump_computeeuv_inputs.pro` (IDL-side input dump)
- `tests/CompareComputeEUVInputs.py` (field-by-field Python vs IDL input comparison)

Additional internal utilities:
- `tests/CompareRenderInputs.py`
- `tests/CompareRenderOutputs.py`
- `tests/CompareModelSources.py`
- `tests/RegressionRoundTripSavH5.py`
- `tests/BuildH5FromSavGroundTruth.py`
