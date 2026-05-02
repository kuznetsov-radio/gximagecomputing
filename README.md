# gximagecomputing / pyGXrender

This code computes 2D maps of the solar microwave (gyroresonance and free-free) and EUV (spectral lines) emission, using models of active regions created by the [**GX Simulator**](https://github.com/Gelu-Nita/GX_SIMULATOR/) (requires SolarSoft GX_Simulator package).

[![Build And Publish Wheels](https://github.com/kuznetsov-radio/gximagecomputing/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/kuznetsov-radio/gximagecomputing/actions/workflows/build_wheels.yml)

----

## PyPI Package Name

PyPI distribution name: `pyGXrender`

Python import package:

```python
import gxrender
```

Install from PyPI:

```bash
pip install pyGXrender
```

----

## Quick Start

See the example files:

- `./examples/idl/RenderExampleMW.pro`
- `./examples/idl/RenderExampleEUV.pro`
- `./examples/python/cli/RenderExampleMW.py`

Examples folder layout:

```text
examples/
├── idl/
│   ├── RenderExampleMW.pro
│   ├── RenderExampleEUV.pro
│   ├── InterpolateEBTELexample.pro
│   └── compile_local_idl
└── python/
    ├── cli/
    │   ├── RenderExampleMW.py
    │   └── RenderExampleEUV.py
    └── sdk/
        ├── sdk_render_mw.py
        └── sdk_render_euv.py
```

**Note:** Sample GX Simulator model and EBTEL data are not included.

----

## Python Environment Notes

### Imports (`PYTHONPATH`)

If running scripts directly from the repository checkout (without installing the package), add `src` to `PYTHONPATH`:

```bash
PYTHONPATH=src python examples/python/cli/RenderExampleMW.py --help
```

If installed with pip, this is usually not needed:

```bash
pip install .
```

### Windows editable installs (`pip install -e .`)

By default, local editable installs on Windows skip native extension compilation
and use prebuilt runtime libraries (for example from `binaries/RenderGRFF_64.dll`).
This avoids requiring Microsoft C++ Build Tools for routine development/testing.

If you explicitly want to compile native code locally, set:

```bash
set PYGXRENDER_BUILD_NATIVE=1
pip install -e .
```

Release wheel builds still force native compilation in CI.

### Writable Config/Cache Directories (SunPy/Matplotlib)

Some environments have non-writable default user config/cache folders. In that case, use writable overrides:

```bash
SUNPY_CONFIGDIR=/tmp/sunpy_cfg MPLCONFIGDIR=/tmp/mpl_cfg python examples/python/cli/RenderExampleMW.py --help
```

If needed, create those folders first:

```bash
mkdir -p /tmp/sunpy_cfg /tmp/mpl_cfg
```

This avoids runtime errors such as:
- `Could not write to SUNPY_CONFIGDIR=...`
- Matplotlib/fontconfig cache permission warnings

## Generated Docs (Sphinx, Optional)

A lightweight Sphinx docs pipeline is available for generated API/reference
pages (SDK, CLI workflow modules, viewer module).

Build locally:

```bash
pip install -r docs/requirements.txt
make docs-html
```

Open:

- `docs/_build/html/index.html`

### EBTEL Table Path (`GXIMAGECOMPUTING_EBTEL_PATH`)

The shared Python workflow/API now requires the EBTEL choice to be explicit:
- pass a real `.sav` path to use DEM/DDM tables
- pass `""` to disable DEM/DDM tables and use the native-library fallback path

The example CLI frontends still honor `GXIMAGECOMPUTING_EBTEL_PATH` as a convenience fallback, but they emit a warning when they do so.

If you want to use that example-layer convenience, set one environment variable once per shell session:

```bash
export GXIMAGECOMPUTING_EBTEL_PATH="$SSW/packages/gx_simulator/euv/ebtel/ebtel_ss.sav"
```

If your SolarSoft installation does not define `$SSW`, use an absolute path:

```bash
export GXIMAGECOMPUTING_EBTEL_PATH="/full/path/to/ssw/packages/gx_simulator/euv/ebtel/ebtel_ss.sav"
```

Then run an example CLI:

```bash
python examples/python/cli/RenderExampleMW.py --model-path /path/to/your.chr.sav --model-format auto
```

### MW Rendering: CLI and Programmatic Usage

`gxrender-mw` is an example-oriented CLI entrypoint. If you omit science inputs such as the frequency list, pixel scale, plasma scalars, or `SHtable`, it will fill the repository's demonstration defaults and emit warnings describing each assumption.

The shared workflow function `gxrender.workflows.render_mw.run(...)` and the SDK entrypoint `gxrender.render_mw_maps(...)` are stricter: they expect those science inputs to be explicit.

If you installed the package (`pip install .` or `pip install -e .`), use the
installed CLI:

```bash
gxrender-mw \
  --model-path /path/to/your.chr.h5 \
  --ebtel-path "" \
  --dx 2.0 --dy 2.0 \
  --frequencies-ghz 5.8 6.2 6.6 7.0 \
  --tbase 1.0e6 --nbase 1.0e8 --q0 0.0217 --a 0.3 --b 2.7 \
  --observer solar-orbiter \
  --output-dir /tmp \
  --output-format h5
```

Optional EBTEL-enabled run:

```bash
gxrender-mw \
  --model-path /path/to/your.chr.h5 \
  --ebtel-path /full/path/to/ebtel_ss.sav \
  --output-dir /tmp \
  --output-format h5
```

You can also call the same workflow from Python:

```python
from argparse import Namespace
from gxrender.workflows.render_mw import run

args = Namespace(
    model_path="/path/to/your.chr.h5",
    model_format="auto",
    ebtel_path="",
    output_dir="/tmp",
    output_name=None,
    output_format="h5",
    omp_threads=None,
    xc=None,
    yc=None,
    observer="earth",
    dx=2.0,
    dy=2.0,
    pixel_scale_arcsec=None,
    nx=None,
    ny=None,
    xrange=None,
    yrange=None,
    frequencies_ghz=[5.8, 6.2, 6.6, 7.0],
    tbase=1.0e6,
    nbase=1.0e8,
    q0=0.0217,
    a=0.3,
    b=2.7,
    corona_mode=0,
    shtable=None,
    shtable_path=None,
    force_isothermal=False,
    interpol_b=False,
    analytical_nt=False,
)
run(args)
```

### Professional SDK Usage (Programmatic, No CLI/argparse Coupling)

For application integration, prefer the SDK layer in `gxrender.sdk` (also re-exported at package root).
This avoids `argparse.Namespace`-style calls and provides typed option objects for MW and EUV rendering.

The SDK does not inject science defaults anymore. Callers are expected to provide explicit plasma scalars, MW frequencies, map pixel scale, and an explicit EBTEL choice (`path` or `""`). The example scripts under `examples/python/` remain the place where demonstration defaults are supplied with warnings.

Available SDK entry points:

- `gxrender.render_mw_maps(...)`
- `gxrender.render_euv_maps(...)`
- `gxrender.MWRenderOptions`
- `gxrender.EUVRenderOptions`
- `gxrender.MapGeometry`
- `gxrender.ObserverOverrides`
- `gxrender.CoronalPlasmaParameters`

MW SDK example:

```python
from gxrender import (
    CoronalPlasmaParameters,
    MapGeometry,
    MWRenderOptions,
    ObserverOverrides,
    render_mw_maps,
)

result = render_mw_maps(
    MWRenderOptions(
        model_path="/path/to/model.chr.sav",
        model_format="sav",
        ebtel_path="/path/to/ebtel.sav",
        output_dir="/tmp/gxrender",
        output_format="h5",
        freqlist_ghz=[5.8, 6.4, 7.2, 8.0, 10.0, 12.0],
        plasma=CoronalPlasmaParameters(
            tbase=1.5e6,
            nbase=2.5e8,
            q0=0.03,
            a=0.4,
            b=2.2,
            mode=0,
            shtable=[
                [1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4],
                [1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4],
                [1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4],
                [1.1, 1.1, 1.1, 1.21, 1.32, 1.43, 1.54],
                [1.2, 1.2, 1.2, 1.32, 1.44, 1.56, 1.68],
                [1.3, 1.3, 1.3, 1.43, 1.56, 1.69, 1.82],
                [1.4, 1.4, 1.4, 1.54, 1.68, 1.82, 0.1],
            ],
        ),
        geometry=MapGeometry(dx=2.0, dy=2.0),
        observer=ObserverOverrides(dsun_cm=None, lonc_deg=None, b0sun_deg=None),
        save_outputs=False,   # in-memory workflow (no files written)
        write_preview=False,  # ignored when save_outputs=False, but explicit is clearer
        verbose=False,  # no CLI-style printing
    )
)

print(result.outputs.h5_path)      # None (save_outputs=False)
print(result.freqlist_ghz[:3])
ti_cube = result.ti                # [ny, nx, nf]
tv_cube = result.tv                # [ny, nx, nf]
print(result.geometry.nx, result.geometry.ny)
```

EUV SDK example:

```python
from gxrender import CoronalPlasmaParameters, EUVRenderOptions, MapGeometry, ObserverOverrides, render_euv_maps

result = render_euv_maps(
    EUVRenderOptions(
        model_path="/path/to/model.chr.sav",
        model_format="sav",
        ebtel_path="/path/to/ebtel.sav",
        response_sav="/path/to/resp_aia_20251126T153431.sav",
        output_dir="/tmp/gxrender",
        geometry=MapGeometry(dx=2.0, dy=2.0),
        plasma=CoronalPlasmaParameters(
            tbase=1.5e6,
            nbase=2.5e8,
            q0=0.03,
            a=0.4,
            b=2.2,
        ),
        observer=ObserverOverrides(
            dsun_cm=14763359700479.328,
            lonc_deg=-17.0574058213,
            b0sun_deg=1.4406505929155138,
        ),
        save_outputs=False,
        write_preview=False,
        verbose=False,
    )
)

print(result.outputs.h5_path)     # None (in-memory render)
print(result.response.channels)
flux_corona = result.flux_corona  # [ny, nx, nch]
flux_tr = result.flux_tr          # [ny, nx, nch]
print(result.outputs.save_outputs, result.outputs.write_preview)
```

Notes:

- The SDK reuses the same rendering engines as the CLI workflows, so CLI and SDK behavior stay aligned.
- CLI entry points remain useful for quick tests and demonstrations.
- The SDK returns typed dataclasses (`MWRenderResult`, `EUVRenderResult`) for a stronger contract than raw dictionaries.
- Set `save_outputs=False` for fully in-memory rendering; set `write_preview=False` to skip preview PNG generation.
- `write_preview` is only used when `save_outputs=True`.
- MW CLI users can override the default science parameters directly with:
  - `--frequencies-ghz`
  - `--tbase`
  - `--nbase`
  - `--q0`
  - `--a`
  - `--b`
  - `--selective-heating`
  - `--shtable-path`
  - `--corona-mode`
  - `--force-isothermal`
  - `--interpol-b`
  - `--analytical-nt`
- `CoronalPlasmaParameters` is shared by MW and EUV and includes the selective-heating connectivity table (`shtable`), so the same plasma/heating configuration can be used in both render paths.
- By default no selective-heating table is applied. To enable selective heating without providing a custom table, set `--selective-heating` in the CLI or `selective_heating=True` in `CoronalPlasmaParameters`; the standard 7x7 table will then be used with a warning.

### Observer Metadata Overrides (MW and EUV, Python CLI)

For parity/debugging workflows, both render CLIs support explicit overrides for
observer/model metadata before calling the native DLL/shared library:

- `--observer` (resolved with SunPy at the true model observation time)
- `--dsun-cm`
- `--lonc-deg`
- `--b0sun-deg`

Observer resolution priority is:

1. CLI observer name (`--observer`)
2. Saved observer state in the model metadata
   - current `pyAMPP` files: `observer/ephemeris/*` first, then `observer/pb0r/*`
   - older headers: `HGLN/HGLT/DSUN` or `CRLN/CRLT/DSUN`
3. Observer name stored in the model metadata
4. CLI triad overrides (`--dsun-cm`, `--lonc-deg`, `--b0sun-deg`) applied on top of the restored/default observer
5. Earth fallback

The resolved geometry is applied before automatic center/FOV inference, so it also
affect default `xc/yc` and FOV calculations unless you pass explicit map
geometry (`--xc`, `--yc`, `--dx`, `--dy`, `--nx`, `--ny`, etc.).

MW example:

```bash
gxrender-mw \
  --model-path /path/to/your.chr.sav \
  --ebtel-path /path/to/ebtel.sav \
  --observer stereo-a
```

EUV example:

```bash
python examples/python/cli/RenderExampleEUV.py \
  --model-path /path/to/your.chr.sav \
  --model-format sav \
  --ebtel-path /path/to/ebtel.sav \
  --response-sav /path/to/resp_aia_20251126T153431.sav \
  --observer "solar orbiter"
```

### Repository Test Fixtures

Large review/test fixtures are distributed separately from the code repository in:

- `https://github.com/suncast-org/pyGXrender-test-data`

Recommended layout:

```text
@SUNCAST-ORG/
  gximagecomputing/
  pyGXrender-test-data/
```

Install the default dataset with:

```bash
cd ../pyGXrender-test-data
scripts/install_dataset.sh
```

The Python shell examples and data-dependent test utilities will auto-detect fixtures from:

```text
../pyGXrender-test-data/raw
```

or from `GXRENDER_TEST_DATA_ROOT` if you prefer a different location.

The default installer populates:

- `raw/models/`
- `raw/responses/`
- `raw/ebtel/`

### Fixture Provenance and Regeneration

The published fixture dataset is reproducible. The model, response, and EBTEL inputs can be regenerated independently.

#### EBTEL tables

The packaged EBTEL tables are also available from the upstream GX Simulator / SolarSoft distribution. In a standard SSW installation they are typically located under:

```text
$SSW/packages/gx_simulator/euv/ebtel/
```

The external fixture repository provides a packaged copy for reproducible testing, but you may also point the examples and tests at your own local SSW/GX Simulator installation.

#### SAV CHR model fixture

The SAV fixture can be regenerated with the original IDL `gx_fov2box` command stored inside the file metadata:

```idl
gx_fov2box, '26-Nov-25 15:47:52', CENTER_ARCSEC=[ -280, -230], DX_KM= 1400, EUV= 1, OUT_DIR='/Users/gelu/Library/CloudStorage/Dropbox/@Projects/sim4fasr/gx_models', SIZE_PIX=[ 150, 100, 100], TMP_DIR='/Users/gelu/Library/CloudStorage/Dropbox/@Projects/sim4fasr/jsoc_cache', UV= 1, CEA= 1
```

This requires an SSW/IDL environment with GX Simulator installed.

#### HDF CHR model fixture

The HDF fixture can be regenerated from the Python `gx-fov2box` command stored in the HDF metadata:

```bash
gx-fov2box --time 2025-11-26T15:47:52 --coords -280.0 -230.0 --hpc --cea --box-dims 150 100 100 --dx-km 1400.000000 --pad-frac 0.1000 --data-dir /Users/gelu/Library/CloudStorage/Dropbox/@Projects/sim4fasr/jsoc_cache --gxmodel-dir /Users/gelu/Library/CloudStorage/Dropbox/@Projects/sim4fasr/gx_models --euv --uv --save-potential --save-bounds --save-nas --save-gen --save-chr --observer-name earth --stop-after chr
```

This requires a working `pyAMPP` / `gx-fov2box` installation and access to the input caches referenced by the command.

#### IDL response files

The response fixtures in `pyGXrender-test-data` were generated in IDL for the test-model epoch using:

- `idlcode/LoadEUVresponse.pro`
- `local/GenerateTestEUVResponses.pro`

The local helper loops over supported instruments and writes date-tagged files such as:

- `resp_aia_20251126T153431.sav`

To regenerate them:

```idl
@/path/to/gximagecomputing/examples/idl/compile_local_idl
.compile '/path/to/gximagecomputing/local/GenerateTestEUVResponses.pro'
GenerateTestEUVResponses
```

This requires an SSW/IDL environment with GX Simulator and the relevant SolarSoft response routines installed. The response generation time is tied to the test-model observation epoch (`2025-11-26T15:34:31`), even though the original `gx_fov2box` request time stored in the model provenance is `2025-11-26T15:47:52`.

### Render Map Viewer GUI (`gxrender-map-view`)

Interactive viewer for rendered MW and EUV map products.

Supported input formats:

- Python-rendered HDF5 map containers (`.h5`, `.hdf5`)
  - MW schema: `maps/data` + `maps/freqlist_ghz`
  - EUV schema: `maps/data` + `maps/channel_ids` (+ optional `maps/component_ids`)
- IDL-rendered map containers (`.sav`, `.xdr`)
  - Combined `map` containers
  - EUV `mapcorona` / `maptr` style containers

CLI usage:

```bash
gxrender-map-view /path/to/rendered_maps.h5
```

Optional initial index (frequency index for MW, channel index for EUV):

```bash
gxrender-map-view /path/to/rendered_maps.h5 --start-index 0
```

Optional solar grid spacing (defaults to 10 degrees; use `0` to disable):

```bash
gxrender-map-view /path/to/rendered_maps.h5 --grid-deg 5
```

Examples:

MW HDF5 output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5_py_mw_maps.h5
```

EUV HDF5 output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.GEN.CHR.h5_py_euv_maps.h5
```

IDL MW SAV output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.sav_idl_mw_maps.sav
```

IDL EUV SAV output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/hmi.M_720s.20201126_195831.E18S19CR.CEA.NAS.CHR.sav_idl_euv_maps.sav
```

Viewer behavior:

- Displays two synchronized panels (left/right components)
  - MW: `TI` and `TV`
  - EUV: `GX (TR)` and `GX (Corona)` (normalized to this order when possible)
- For Python-rendered HDF5, reads the stored SunPy-ready WCS header from `metadata/wcs_header`
- Performs no observer or geometry reconstruction for HDF5 inputs
- Uses map-appropriate units
  - MW: `K`
  - EUV: `DN s^-1 pix^-1`
- Replaces `NaN/Inf` pixels with `0` for robust display
- Draws a solar coordinate grid by default (`--grid-deg 10`)
- Provides per-panel controls:
  - intensity range slider
  - log scaling toggle
- Provides a shared axis slider:
  - MW: frequency index
  - EUV: channel index

Notes:

- For EUV HDF5 files, the slider moves across channels (e.g. `A94`, `A131`, ...)
- For IDL SAV files, the viewer auto-detects MW vs EUV from map IDs
- If your environment has restrictive config/cache permissions, set:

```bash
SUNPY_CONFIGDIR=/tmp/sunpy_cfg MPLCONFIGDIR=/tmp/mpl_cfg gxrender-map-view /path/to/rendered_maps.h5
```

----

## Python Data Branches (Non-Interfering)

The Python API now treats CHR inputs as two explicit branches that both normalize to the same internal `ChromoModel` representation consumed by the rendering library:

- **IDL branch**: `load_model_sav(...)` for GX Simulator `.sav` CHR models.
- **pyAMPP branch**: `load_model_hdf(...)` for current pyAMPP `.h5` CHR models (`/chromo` group).

Both branches are converted into one internal data layout before calling the native renderer, so loader-specific format changes do not leak into rendering logic.

### Loader Observer Contract

The Python and IDL model loaders now follow the same observer-precedence contract:

- explicit scalar overrides win:
  - `DSun`
  - `lonC`
  - `b0Sun`
- otherwise, if `recompute_observer_ephemeris` is requested, the loader recomputes observer ephemeris from time using the requested observer name, or the saved observer name when it is a standard resolvable observer
- otherwise, the loader uses observer geometry recovered from the file itself
- only if the file does not contain usable observer geometry does the loader fall back to implicit Earth defaults

This means saved non-Earth HDF5/SAV products are loaded as-saved by default. Recompute mode is an explicit opt-in reinterpretation step.

For `pyAMPP` HDF5 inputs, the file-first observer geometry follows the saved `pyAMPP` convention directly:

- `observer/pb0r/b0_deg` and `observer/pb0r/l0_deg` are treated as the saved `B0/L0` observer angles
- these are true heliographic Stonyhurst observer angles
- in particular, `observer/pb0r/l0_deg` is the saved Stonyhurst `L0`, not `CRLN_OBS` and not the model `lonC`
- `observer/pb0r/rsun_arcsec` is the saved apparent solar radius for that observer geometry

So `gxrender` does not define a new saved-file observer convention here; it follows the one already written by `pyAMPP`.

On the IDL HDF path, this is handled by `ConvertToGX.pro` preserving the saved HDF5 `/observer` subtree as a separate nested `box.observer` structure, adapting dynamically to whatever keys are present. The original `box.index` header is treated as the saved birth-certificate header and is not rewritten from `/observer/*`.

If a caller needs both the strict DLL-ready model structure and the saved
observer metadata without reading the input model twice:

- IDL: `LoadGXmodel, modelfile, observer_struct=observer, index_struct=index`
- Python HDF5: `load_model_hdf_with_observer(...)`
- Python SAV: `load_model_sav_with_observer(...)`

These return the normal DLL-ready model plus the saved observer metadata group
as stored in the input file. This keeps saved LOS state and saved FOV metadata
available for later image-geometry decisions without changing the DLL model or
simbox conventions.

For Python, the returned tuple is:

- `model`: the unchanged DLL-ready model structure
- `model_dt`: the NumPy dtype describing that structure
- `header`: resolved loader metadata after file-first defaults, recompute, and/or explicit overrides
- `observer`: the saved observer metadata group as stored in the file

For the IDL render examples, non-Earth LOS and automatic FOV handling now run
through the shared geometry helper suite in `idlcode/GXObserverGeometry.pro`:

- `GXResolveObserverGeometry`
- `GXComputeInscribingFOV`
- `GXResolveSimboxFromObserverAndModel`

This keeps the DLL-facing `model` and `simbox` conventions unchanged while
making the IDL path follow the same saved-observer and saved-FOV logic as the
Python renderer.

By default, the IDL render examples use the saved `observer/fov` rectangle when
it is present and no explicit observer or view overrides were requested. To
force a fresh inscribing-FOV computation instead, pass `/AUTO_FOV` to the IDL
render example. `USE_SAVED_FOV` remains available as a backward-compatible
alias, but the preferred interface is:

- default: use saved `observer/fov` when present
- `/AUTO_FOV`: recompute the observer-aligned inscribing FOV from geometry

Parity expectations for this contract are:

- `default` file-first loading is expected to be numerically equivalent between the Python and IDL loaders for the same saved input model
- explicit scalar overrides (`DSun`, `lonC`, `b0Sun`) are expected to be numerically equivalent between the Python and IDL loaders
- `recompute_observer_ephemeris` is not expected to be numerically identical between Python and IDL, because each side uses its own ephemeris engine

In other words, recompute mode is a deliberate observer reinterpretation feature, not a cross-language bitwise parity mode. For reproducible parity checks, use the saved file observer state or explicit scalar overrides.

### Native SAV -> HDF5 Conversion (No pyAMPP Dependency)

Use the built-in CLI to convert a GX CHR `.sav` model into canonical HDF5:

```bash
gx-sav2h5 \
  --sav-path /path/to/input.NAS.CHR.sav \
  --out-h5 /path/to/output.NAS.CHR.h5
```

Optional: seed from an existing HDF5 template while still rewriting model groups:

```bash
gx-sav2h5 \
  --sav-path /path/to/input.NAS.CHR.sav \
  --out-h5 /path/to/output.NAS.CHR.h5 \
  --template-h5 /path/to/template.h5
```

### Internal Validation Workflows

Repository-internal parity/regression procedures (including IDL/Python parity and
comparison scripts under `tests/`) are documented in `tests/README.md`.

----

## Building Native Library (Linux/macOS)

The `source/makefile` supports platform-aware builds and copies outputs into `./binaries`.

### Build

```bash
cd source
make
```

### Outputs

- Linux: `binaries/RenderGRFF.so`
- macOS arm64: `binaries/RenderGRFF_arm64.so`
- macOS x86_64: `binaries/RenderGRFF_x86_64.so`

### macOS prerequisites

Install OpenMP runtime (Homebrew):

```bash
brew install libomp
```

If Homebrew is in a non-default prefix, set include/link flags explicitly:

```bash
make CPPFLAGS='-I/opt/homebrew/opt/libomp/include' LDFLAGS='-L/opt/homebrew/opt/libomp/lib'
```

### Binary Wheel Releases

For maintainers: release process and exact publish commands are documented in `RELEASING.md`.
CI workflow: `.github/workflows/build_wheels.yml`

----

## Microwave Emission Maps

### Step 1: Load the GX Simulator Model

```idl
model = LoadGXmodel(modelfile [, newTime=newTime])
```
- `modelfile`: GX Simulator model file name (must include field line info and chromospheric part).
- `newTime`: Optional date/time (accepted by `anytim()`), rotates the model to new date/time if specified.

---

### Step 2: Load the EBTEL Tables

```idl
ebtel = LoadEBTEL(ebtelfile [, DEM=DEM, DDM=DDM])
```
- `ebtelfile`: GX Simulator sav file with EBTEL table(s) (DEM and/or DDM).
- If `ebtelfile=''`: DEM, DDM and coronal heating model are **not used** (coronal plasma described by constant temperature and barometric height profile of density).
- `/DEM` and `/DDM` keywords: Use when both DEM and DDM tables are present.
- If only one table exists, keywords are ignored.

---

### Step 3: Define Map Size, Position and Frequencies

```idl
simbox = MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, freqlist [, rot=rot, /parallel, /exact, Nthreads=Nthreads])
```

- `xc, yc`: Center (helioprojective x, y) in arcseconds.
- `dx, dy`: Output resolution, arcseconds.
- `Nx, Ny`: Map size (pixels).
- `freqlist`: 1D array of emission frequencies (GHz).
- `rot`: Optional rotation angle (degrees, counterclockwise). The center of the resulting map is still at (xc, yc), and the x and y coordinates correspond to the rotated coordinate system
- `/parallel`: Render with parallel projection (all lines of sight are parallel to each other) (default: perspective, all lines of sight intersect at the observer's location).
- `/exact`: Use with `/parallel`, controls conversion to kilometers. If not set (default), the conversion from arcseconds to kilometers in the parallel projection is performed using the distance from the observer to the center of the Sun. If set, the conversion is performed using the actual distance from the observer to the considered active region.
- `Nthreads`: Number of processor threads (≤ available processors). Default: a system-defined value (typically, the number of available processors).

---

### Step 4: Define Coronal Plasma Parameters

```idl
coronaparms = DefineCoronaParams(Tbase, nbase, Q0, a, b [, /force_isothermal, /analyticalNT])
```
- `Tbase`: Plasma temperature (K).
- `nbase`: Base plasma density (cm^{-3}) at the bottom of the simulation box.
- `Q0, a, b`: Coronal heating model (applies to closed field lines), where heating rate:
  ```
  Q = Q0*(B/B0)^a / (L/L0)^b
  ```
- `/force_isothermal`: Ignore multi-thermal formulae given in the paper of Fleishman, Kuznetsov & Landi (2021), use the moments of the DEM or DDM distribution (if both DEM and DDM are provided, the DDM moments are used). This option improves the computation speed greatly, although the results become less accurate.
- `/analyticalNT`: Use analytical formula for voxels with heating parameters outside EBTEL table bounds.

`Tbase` and `nbase` are used to find the plasma parameters in the voxels associated with open field lines, or, if the keyword `/analyticalNT` is not set, the heating parameters (Q, L) in closed field lines which are beyond the boundaries of the EBTEL table. In such voxels, the plasma temperature is set to `Tbase`, and the plasma density is computed using `nbase`, `Tbase`, and the barometric formula.

---

### Step 5: Prepare Output Memory Structure

```idl
outspace = ReserveOutputSpace(simbox)
```

---

### Step 6: (Optional) Selective Heating Table

Prepare `SHtable` (2D double array) to define selective heating for coronal magnetic field lines.

---

### Main Microwave Computation

Call the main executable module (RenderGRFF_32.dll, RenderGRFF_64.dll, or RenderGRFF.so) via `call_external`:

```idl
r = call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace [, SHtable])
```

- `libname`: Name of executable library
- Remaining arguments: Structures from above steps

#### Output Structure

- `outspace.TI` & `outspace.TV`: Brightness temperatures for Stokes parameters I and V (in K). Each field is a 3D array with Nx * Ny * Nf elements, where Nx and Ny are the x and y sizes of the computed maps, and Nf is the number of the emission frequencies.
- To convert:
    ```idl
    ConvertToMaps, outspace, simbox, model, mapI, mapV [, /flux]
    ```
    - `mapI`, `mapV`: SolarSoft multi-frequency map objects
    - If `/flux` is specified, unit changes to sfu/pix
- Additional fields: `outspace.flagsAll`, `outspace.flagsCorona` (see [Computation Statistics](#computation-statistics)).

---

#### Example Usage

See: `examples/idl/RenderExampleMW.pro` (sample data not included).

----

## EUV Emission Maps

Steps are similar, but with differences noted below.

### Step 1: Load the GX Simulator Model

```idl
model = LoadGXmodel(modelfile [, newTime=newTime])
```
*(see microwave section above for description)*

---

### Step 2: Load the EBTEL Tables

```idl
ebtel = LoadEBTEL(ebtelfile)
```

- Includes EBTEL tables for DEM (corona and transition region).
- If `ebtelfile=''`: DEM and heating model **not used**.

The keyword `/DDM` should not be used, because the EUV emission depends on the DEM only.

---

### Step 3: Load Instrumental Response Function

```idl
response = LoadEUVresponse(model.obstime [, instrument, evenorm=evenorm, chiantifix=chiantifix])
```
- `model.obstime`: Observation time from `LoadGXmodel`
- `instrument`: Choose from `'AIA'`, `'AIA2'`, `'TRACE'`, `'SXT'`, `'SOLO-FSI'`, `'SOLO-HRI'`, `'STEREO-A'`, `'STEREO-B'` (default `'AIA'`).
- `evenorm`, `chiantifix`: AIA parameters, default=1 (see SolarSoft `aia_get_response.pro`).

For backward compatibility, `LoadEUVresponse` also accepts a full model structure and will read its `OBSTIME` field, but the preferred interface is to pass the time directly.

---

### Step 4: Define EUV Map Size and Position

```idl
simbox = MakeSimulationBoxEUV(xc, yc, dx, dy, Nx, Ny [, /parallel, /exact, Nthreads=Nthreads])
```

- Channels: All specified by instrument's response table (cannot select individual channels).
- Emission computed as observed from Earth's distance. Thus for Solar Orbiter and STEREO the map position and pixel size should be corrected accordingly

---

### Step 5: Define Coronal Plasma Parameters

```idl
coronaparms = DefineCoronaParams(Tbase, nbase, Q0, a, b [, /analyticalNT])
```
- `/force_isothermal` has **no effect** for EUV emission.
- Other parameters, see microwave emission above.

---

### Step 6: Prepare Output Memory Structure

```idl
outspace = ReserveOutputSpaceEUV(simbox, response)
```
- `simbox` from above, `response` from LoadEUVresponse

---

### Step 7: (Optional) Selective Heating Table

Prepare `SHtable` as described above.

---

### Main EUV Computation

Call the main executable module as:

```idl
r = call_external(libname, 'ComputeEUV', model, ebtel, response, simbox, coronaparms, outspace [, SHtable])
```

#### Output Structure

- `outspace.fluxCorona`, `outspace.fluxTR`: Computed EUV fluxes (units: DN s^{-1} pix^{-1})
- To convert:
    ```idl
    ConvertToMapsEUV, outspace, simbox, model, response, mapCorona, mapTR
    ```
    - `mapCorona`, `mapTR`: SolarSoft multi-channel maps

- Flags information: `outspace.flagsAll`, `outspace.flagsCorona` (see [Computation Statistics](#computation-statistics))

#### Example Usage

See: `examples/idl/RenderExampleEUV.pro` (sample data not included).

----

## Selective Heating Table (`SHtable`)

Both microwave and EUV emission computations may utilize the optional selective heating table (`SHtable`):

- A 2D array (double precision), with 7 * 7 elements:
    ```
    [number of closed field lines, number of simulation epochs]
    ```
- Default value: `1.0` for all of the elements. Each element of that table represents the factor applied to the heating rate Q for the field lines connecting specific regions at the photosphere; see the 'Selective Heating Mask' panel in GX Simulator. The SHtable table is supposed to be symmetric, i.e., `SHtable[j, i]=SHtable[i, j]`; asymmetric tables are accepted but the result will likely have no sense.
- Used in both `ComputeMW` and `ComputeEUV` calls when provided.

----

## <a name="computation-statistics"></a>Computation Statistics: Output Flags

The output structure contains fields for computation statistics:

### `flagsAll` (length=6):

| Index | Meaning                                                              |
|-------|-----------------------------------------------------------------------|
| 0     | Total number of voxels crossed by lines-of-sight                      |
| 1     | Number of voxels in chromospheric part of model (crossed)             |
| 2     | Number of voxels (crossed by the LOS) associated with closed field lines (known loop length `L` and average magnetic field `B_avg`). `flagsAll[2]=flagsAll[3]+flagsAll[4]+flagsAll[5]` |
| 3     | Voxels (crossed and closed field lines) with EBTEL table hits (both `L` and `Q` within table) |
| 4     | Voxels (crossed and closed field lines) missing EBTEL table due to loop length (`L` is beyond the table) |
| 5     | Voxels (crossed and closed field lines) missing EBTEL table due to heating rate (`Q` out of bounds) |

### `flagsCorona` (length=6):

Similar to `flagsAll`, but refers only to the *coronal* part of the model.

- `flagsCorona[0]`: Total number of voxels in the coronal part crossed by lines-of-sight
- `flagsCorona[1]`: always zero

----

## References

For detailed theory and formulae, see the relevant publications, especially:
- [Fleishman, Kuznetsov & Landi (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...914...52F/abstract)

For questions or issues, please open a GitHub issue or contact the author.

----
