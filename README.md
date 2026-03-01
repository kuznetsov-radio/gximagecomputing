# gximagecomputing

This code computes 2D maps of the solar microwave (gyroresonance and free-free) and EUV (spectral lines) emission, using models of active regions created by the [**GX Simulator**](https://github.com/Gelu-Nita/GX_SIMULATOR/) (requires SolarSoft GX_Simulator package).

[![Build And Publish Wheels](https://github.com/kuznetsov-radio/gximagecomputing/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/kuznetsov-radio/gximagecomputing/actions/workflows/build_wheels.yml)

----

## PyPI Package Name

PyPI distribution name: `pyGXrender`

Python import package remains:

```python
import gximagecomputing
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

EBTEL is optional. If no EBTEL path is provided, DEM/DDM and heating tables are not used (legacy behavior: isothermal/hydrostatic fallback in the native library).

If you want to use EBTEL tables, pass `--ebtel-path` explicitly or set one environment variable once per shell session:

```bash
export GXIMAGECOMPUTING_EBTEL_PATH="$SSW/packages/gx_simulator/euv/ebtel/ebtel_ss.sav"
```

If your SolarSoft installation does not define `$SSW`, use an absolute path:

```bash
export GXIMAGECOMPUTING_EBTEL_PATH="/full/path/to/ssw/packages/gx_simulator/euv/ebtel/ebtel_ss.sav"
```

Then run:

```bash
python examples/python/cli/RenderExampleMW.py --model-path /path/to/your.chr.sav --model-format auto
```

### MW Rendering: CLI and Programmatic Usage

If you installed the package (`pip install .` or `pip install -e .`), use the
installed CLI:

```bash
gxrender-mw \
  --model-path /path/to/your.chr.h5 \
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
from gximagecomputing.workflows.render_mw import run

args = Namespace(
    model_path="/path/to/your.chr.h5",
    model_format="auto",
    ebtel_path=None,
    output_dir="/tmp",
    output_name=None,
    output_format="h5",
    omp_threads=None,
    xc=None,
    yc=None,
    dx=2.0,
    dy=2.0,
    pixel_scale_arcsec=None,
    nx=None,
    ny=None,
    xrange=None,
    yrange=None,
)
run(args)
```

### Professional SDK Usage (Programmatic, No CLI/argparse Coupling)

For application integration, prefer the SDK layer in `gximagecomputing.sdk` (also re-exported at package root).
This avoids `argparse.Namespace`-style calls and provides typed option objects for MW and EUV rendering.

Available SDK entry points:

- `gximagecomputing.render_mw_maps(...)`
- `gximagecomputing.render_euv_maps(...)`
- `gximagecomputing.MWRenderOptions`
- `gximagecomputing.EUVRenderOptions`
- `gximagecomputing.MapGeometry`
- `gximagecomputing.ObserverOverrides`

MW SDK example:

```python
from gximagecomputing import MapGeometry, MWRenderOptions, ObserverOverrides, render_mw_maps

result = render_mw_maps(
    MWRenderOptions(
        model_path="/path/to/model.chr.sav",
        model_format="sav",
        ebtel_path="/path/to/ebtel.sav",
        output_dir="/tmp/gxrender",
        output_format="h5",
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
from gximagecomputing import EUVRenderOptions, MapGeometry, ObserverOverrides, render_euv_maps

result = render_euv_maps(
    EUVRenderOptions(
        model_path="/path/to/model.chr.sav",
        model_format="sav",
        ebtel_path="/path/to/ebtel.sav",
        response_sav="/path/to/aia_response.sav",
        output_dir="/tmp/gxrender",
        geometry=MapGeometry(dx=2.0, dy=2.0),
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

### Observer Metadata Overrides (MW and EUV, Python CLI)

For parity/debugging workflows, both render CLIs support explicit overrides for
observer/model metadata before calling the native DLL/shared library:

- `--dsun-cm`
- `--lonc-deg`
- `--b0sun-deg`

These overrides are applied before automatic center/FOV inference, so they also
affect default `xc/yc` and FOV calculations unless you pass explicit map
geometry (`--xc`, `--yc`, `--dx`, `--dy`, `--nx`, `--ny`, etc.).

MW example:

```bash
gxrender-mw \
  --model-path /path/to/your.chr.sav \
  --ebtel-path /path/to/ebtel.sav \
  --dsun-cm 14763359700479.328 \
  --lonc-deg -17.0574058213 \
  --b0sun-deg 1.4406505929155138
```

EUV example:

```bash
python examples/python/cli/RenderExampleEUV.py \
  --model-path /path/to/your.chr.sav \
  --model-format sav \
  --ebtel-path /path/to/ebtel.sav \
  --response-sav /path/to/aia_response.sav \
  --dsun-cm 14763359700479.328 \
  --lonc-deg -17.0574058213 \
  --b0sun-deg 1.4406505929155138
```

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

Examples:

MW HDF5 output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_py_mw_maps.h5
```

EUV HDF5 output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_py_euv_maps.h5
```

IDL MW SAV output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_idl_mw_maps.sav
```

IDL EUV SAV output:

```bash
gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_idl_euv_maps_forced_observer.sav
```

Viewer behavior:

- Displays two synchronized panels (left/right components)
  - MW: `TI` and `TV`
  - EUV: `GX (TR)` and `GX (Corona)` (normalized to this order when possible)
- Preserves WCS metadata when available from HDF5 `metadata/index_header`
- Uses map-appropriate units
  - MW: `K`
  - EUV: `DN s^-1 pix^-1`
- Replaces `NaN/Inf` pixels with `0` for robust display
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
response = LoadEUVresponse(model [, instrument, evenorm=evenorm, chiantifix=chiantifix])
```
- `model`: From `LoadGXmodel`
- `instrument`: Choose from `'AIA'`, `'AIA2'`, `'TRACE'`, `'SXT'`, `'SOLO-FSI'`, `'SOLO-HRI'`, `'STEREO-A'`, `'STEREO-B'` (default `'AIA'`).
- `evenorm`, `chiantifix`: AIA parameters, default=1 (see SolarSoft `aia_get_response.pro`).

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
