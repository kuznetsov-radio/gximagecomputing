# Tests / Internal Validation Workflows

This document is for repository-internal regression and parity workflows.
These procedures are useful for development/debugging and are not part of the
public user-facing quick-start documentation.

## Contents

- External test-data installation
- SAV↔H5 model parity checks
- MW/EUV renderexample wrappers on Unix, macOS, and Windows
- One-command Python-vs-IDL renderexample parity benchmark driver
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

These wrappers are intentionally thin launchers around the Python workflow CLIs.
They resolve repository test fixtures, set isolated SunPy/Matplotlib cache
directories, and then call `src/gxrender/workflows/render_euv.py` or
`src/gxrender/workflows/render_mw.py`.

Show the full option surface with:

```bash
bash scripts/unix/renderexampleeuv_test.sh --help
bash scripts/unix/renderexamplemw_test.sh --help
```

```bat
scripts\win\renderexampleeuv_test.bat --help
scripts\win\renderexamplemw_test.bat --help
```

Common options shared by the Unix and Windows wrappers:

- `--model-path PATH`: explicit H5/SAV model input
- `--model-name NAME`: named model fixture from the external test-data set
- `--ebtel PATH`: explicit EBTEL table
- `--ebtel-name NAME`: named EBTEL fixture, default `ebtel.sav`
- `--output-dir PATH` and `--output-name NAME`: output location/name
- `--observer NAME`: observer override such as `earth`, `stereo-a`, `stereo-b`, or `solar orbiter`
- `--auto-fov`: recompute the observer-aligned inscribing FOV
- `--use-saved-fov`: force saved-FOV preference
- `--xc`, `--yc`, `--dx`, `--dy`, `--pixel-scale-arcsec`, `--nx`, `--ny`, `--xrange`, `--yrange`: explicit map geometry
- `--dsun-cm`, `--lonc-deg`, `--b0sun-deg`: explicit observer metadata overrides
- `--q0`, `--a`, `--b`: closed-field heating overrides
- `--show-maps`: open `gxrender-map-view` after rendering when available

EUV-only options:

- `--instrument NAME`: EUV instrument override
- `--channels CH0 [CH1 ...]`: explicit EUV channel list
- `--response PATH`: explicit EUV response SAV

MW-only options:

- `--frequencies-ghz F0 [F1 ...]`: explicit MW frequency list
- `--freqlist-ghz F0 [F1 ...]`: alias for `--frequencies-ghz`

Environment variables provide the same defaults when scripting:

- `PYTHON_BIN`: Python interpreter used by the wrapper
- `GXRENDER_TEST_DATA_ROOT`: external fixture root
- `OUTDIR`, `OUTNAME`, `MODEL_PATH`, `MODEL_NAME`, `EBTEL_PATH`, `EBTEL_NAME`
- `OBSERVER`, `AUTO_FOV`, `USE_SAVED_FOV`, `RUNTIME_CACHE_ROOT`
- `MPLCONFIGDIR`, `SUNPY_CONFIGDIR`
- `RESPONSE` for Unix EUV; `RESPONSE` or `RESPONSE_SAV` for Windows EUV

Examples:

```bash
bash scripts/unix/renderexampleeuv_test.sh \
  --model-path /path/to/model.NAS.GEN.CHR.h5 \
  --response /path/to/aia_response.sav \
  --channels 171 193 211 \
  --auto-fov \
  --output-dir /tmp/gximagecomputing_euv
```

```bash
bash scripts/unix/renderexamplemw_test.sh \
  --model-path /path/to/model.NAS.GEN.CHR.h5 \
  --frequencies-ghz 5.8 8.0 10.0 \
  --pixel-scale-arcsec 2.0 \
  --nx 128 --ny 128 \
  --output-dir /tmp/gximagecomputing_mw
```

Windows `cmd.exe` equivalents:

```bat
scripts\win\renderexampleeuv_test.bat ^
  --model-path C:\data\model.NAS.GEN.CHR.h5 ^
  --response C:\data\aia_response.sav ^
  --channels 171 193 211 ^
  --auto-fov ^
  --output-dir %TEMP%\gximagecomputing_euv
```

```bat
scripts\win\renderexamplemw_test.bat ^
  --model-path C:\data\model.NAS.GEN.CHR.h5 ^
  --frequencies-ghz 5.8 8.0 10.0 ^
  --pixel-scale-arcsec 2.0 ^
  --nx 128 --ny 128 ^
  --output-dir %TEMP%\gximagecomputing_mw
```

## One-Command RenderExample Python-vs-IDL Benchmark

`scripts/unix/run_renderexample_parity_benchmarks.sh` is the CI-oriented driver
for end-to-end Python-vs-IDL renderexample parity. It runs the Python wrappers,
generates temporary IDL batch files for `RenderExampleMW_test` and/or
`RenderExampleEUV_test`, compares the output maps, and writes comparison JSON
summaries under the output root.

Basic usage:

```bash
scripts/unix/run_renderexample_parity_benchmarks.sh --help
scripts/unix/run_renderexample_parity_benchmarks.sh --mode both
```

Modes:

- `--mode mw`: MW only
- `--mode euv`: EUV only
- `--mode both`: MW and EUV, the default

Important inputs:

- `--python PATH` or `PYTHON_BIN`: Python interpreter. If omitted, the script uses `python3` or `python` from `PATH`.
- `--idl PATH` or `IDL_BIN`: IDL launcher. If omitted, the script uses `sswidl` or `idl` from `PATH`.
- `--model-h5 PATH` or `MODEL_H5_PATH`: H5 model input for the Python workflow.
- `--model-sav PATH` or `MODEL_SAV_PATH`: SAV model input for the IDL workflow.
- `--ebtel PATH` or `GXIMAGECOMPUTING_EBTEL_PATH`: EBTEL table.
- `--response-sav PATH` or `GXIMAGECOMPUTING_EUV_RESPONSE_SAV`: EUV response SAV, required for `euv` and `both`.
- `--outroot PATH` or `OUTROOT`: artifact root. Defaults to `$RUNNER_TEMP/gximagecomputing_renderexample_parity`, `$TMPDIR/gximagecomputing_renderexample_parity`, or `/tmp/gximagecomputing_renderexample_parity`.

Prefer `MODEL_H5_PATH` plus `MODEL_SAV_PATH` in CI. `MODEL_PATH` / `--model-path`
is kept only as a backward-compatible fallback and uses the same file for both
Python and IDL inputs.

Example local run:

```bash
scripts/unix/run_renderexample_parity_benchmarks.sh \
  --mode both \
  --python "$(command -v python3)" \
  --idl "$(command -v sswidl)" \
  --model-h5 /path/to/model.NAS.GEN.CHR.h5 \
  --model-sav /path/to/model.NAS.CHR.sav \
  --ebtel /path/to/ebtel.sav \
  --response-sav /path/to/aia_response.sav \
  --outroot /tmp/gximagecomputing_renderexample_parity
```

Example CI step:

```yaml
- name: RenderExample parity
  run: |
    scripts/unix/run_renderexample_parity_benchmarks.sh \
      --mode both \
      --model-h5 "$MODEL_H5_PATH" \
      --model-sav "$MODEL_SAV_PATH" \
      --ebtel "$GXIMAGECOMPUTING_EBTEL_PATH" \
      --response-sav "$GXIMAGECOMPUTING_EUV_RESPONSE_SAV" \
      --outroot "$RUNNER_TEMP/renderexample-parity"
  env:
    PYTHON_BIN: ${{ env.pythonLocation }}/bin/python
    IDL_BIN: sswidl
```

Artifacts:

- MW Python output: `$OUTROOT/mw/*_py_mw_maps.h5`
- MW IDL output: `$OUTROOT/mw/*_idl_mw_maps.sav`
- MW comparison JSON: `$OUTROOT/mw/comparison_python_vs_idl.json`
- EUV Python output: `$OUTROOT/euv/*_py_euv_maps.h5`
- EUV IDL output: `$OUTROOT/euv/*_idl_euv_maps.sav`
- EUV comparison JSON: `$OUTROOT/euv/comparison_python_vs_idl.json`

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
