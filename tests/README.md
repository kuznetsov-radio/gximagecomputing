# Tests / Internal Validation Workflows

This document is for repository-internal regression and parity workflows.
These procedures are useful for development/debugging and are not part of the
public user-facing quick-start documentation.

## Contents

- External test-data installation
- SAV↔H5 model parity checks
- IDL/Python MW/EUV map parity comparisons
- ComputeEUV pre-DLL input dump/compare workflows

## External Test-Data Installation

The large CHR model fixtures, EUV response bundles, and EBTEL tables are now distributed separately from this repository.

Recommended layout:

```text
@SUNCAST-ORG/
  gximagecomputing/
  pyGXrender-test-data/
```

Clone the data repository next to `gximagecomputing` and install the default fixture set:

```bash
cd ..
git clone https://github.com/suncast-org/pyGXrender-test-data.git
cd pyGXrender-test-data
scripts/install_dataset.sh
```

By default, the Python tests and workflow wrappers will look for extracted fixtures under:

```text
../pyGXrender-test-data/raw
```

You can override that discovery path with:

```bash
export GXRENDER_TEST_DATA_ROOT=/path/to/pyGXrender-test-data/raw
```

## One-Command SAV↔H5 Parity Check

Run strict regression parity (rebuild H5 from SAV, then compare loader outputs field-by-field):

```bash
make parity-roundtrip
```

Defaults:
- `SAV_PATH` and `H5_PATH` are auto-resolved from the external dataset when omitted
- `ATOL=0`, `RTOL=0`

Override example:

```bash
make parity-roundtrip \
  SAV_PATH=/path/to/model.NAS.CHR.sav \
  H5_PATH=/tmp/model.NAS.CHR.h5
```

If your fixture directory is under Dropbox/iCloud and file locking interferes, prefer a temporary output file:

```bash
make parity-roundtrip H5_PATH=/tmp/gximagecomputing_roundtrip_from_sav.h5
```

## Workflow Wrappers By Platform

Run these from repository root (`gximagecomputing/`) after activating your Python environment.

Unix/macOS (or Git Bash):

```bash
bash scripts/unix/renderexampleeuv_test.sh
bash scripts/unix/renderexamplemw_test.sh
```

Windows cmd/PowerShell:

```bat
scripts\win\renderexampleeuv_test.bat
scripts\win\renderexamplemw_test.bat
```

Windows from Git Bash:

```bash
cmd //c scripts\\win\\renderexampleeuv_test.bat
cmd //c scripts\\win\\renderexamplemw_test.bat
```

Default output location:

- Unix/macOS wrappers: `/tmp/gximagecomputing_validation_groundtruth`
- Windows wrappers: `%TEMP%\gximagecomputing_validation_groundtruth`

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
  RESPonsefile='/path/to/resp_aia_20251126T153431.sav', $
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
  --response-sav /path/to/resp_aia_20251126T153431.sav \
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
PYTHONPATH=src python scripts/python/ComparePythonVsIDLEUVMaps.py \
  --python-h5 /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5_py_euv_maps.h5 \
  --idl-sav /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.sav_idl_euv_maps.sav
```

## Compare Scripts

MW/EUV map comparisons:
- `scripts/python/ComparePythonVsIDLMaps.py` (shared implementation; supports `--kind auto|mw|euv`)
- `scripts/python/ComparePythonVsIDLEUVMaps.py` (EUV wrapper with `--kind euv` default)

ComputeEUV pre-DLL input parity:
- `scripts/python/DumpComputeEUVInputs.py` (Python-side input dump)
- `scripts/idl/dump_computeeuv_inputs.pro` (IDL-side input dump)
- `scripts/python/CompareComputeEUVInputs.py` (field-by-field Python vs IDL input comparison)

Additional internal utilities:
- `scripts/python/CompareRenderInputs.py`
- `scripts/python/CompareRenderOutputs.py`
- `scripts/python/CompareModelSources.py`
- `scripts/python/RegressionRoundTripSavH5.py`
- `scripts/python/BuildH5FromSavGroundTruth.py`
