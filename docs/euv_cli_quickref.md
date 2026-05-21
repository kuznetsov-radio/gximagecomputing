Quick Reference: EUV Preview & Colormap CLI
===============================================

Command: gxrender-euv-preview
-----------------------------

**Purpose:** Regenerate EUV preview PNGs from rendered H5 maps with flexible channel selection and log scaling.

**Syntax:**

    gxrender-euv-preview MAP_PATH [OPTIONS]

**Positional Arguments:**

    MAP_PATH              Path to rendered EUV map H5 file

**Optional Arguments:**

    --channel-index INT      Channel index (0-based). Default: 0
    --channel-id LABEL       Channel label (e.g., A171, A284, 195)
    --output PATH            Output PNG path. Default: <stem>_preview_<channel>.png
    --show                   Display with matplotlib
    --no-save                Suppress PNG file writing
    --title TEXT             Override preview title
    --log-scale              Apply log10 scaling for dynamic range enhancement

**Examples:**

    # Generate linear preview for channel A171 (default)
    gxrender-euv-preview map.h5

    # Generate log-scaled preview for channel A284
    gxrender-euv-preview map.h5 --channel-id A284 --log-scale

    # By index
    gxrender-euv-preview map.h5 --channel-index 2

    # Custom output path
    gxrender-euv-preview map.h5 --channel-id A171 --output /tmp/my_preview.png

    # Display interactively without saving
    gxrender-euv-preview map.h5 --channel-id A195 --no-save --show

    # Check available channels
    python -c "import h5py; f = h5py.File('map.h5'); print(f['maps']['channel_ids'][:])"

----

Command: gxrender-euv
---------------------

**Purpose:** Render EUV maps from CHR model with optional preview generation.

**Syntax:**

    gxrender-euv --model-path MODEL [OPTIONS]

**Core Arguments:**

    --model-path PATH           Path to CHR model (.h5 or .sav)
    --model-format {h5,sav,auto}
                                Format detection. Default: auto
    --ebtel-path PATH           EBTEL table (.sav) or "" to disable
    --output-dir PATH           Output directory. Default: ./
    --output-name NAME          H5 output name. Default: <model>_py_euv_maps.h5

**Geometry Arguments:**

    --dx ARCSEC/PIX             Pixel scale (arcsec/pix). Default: 2.0
    --dy ARCSEC/PIX             Pixel scale (arcsec/pix). Default: 2.0
    --pixel-scale-arcsec FLOAT  Set both dx and dy
    --nx INT                    Force image width (pixels)
    --ny INT                    Force image height (pixels)
    --xc ARCSEC                 Center X (default: model metadata)
    --yc ARCSEC                 Center Y (default: model metadata)
    --xrange X1 X2              X bounds (arcsec)
    --yrange Y1 Y2              Y bounds (arcsec)

**Observer Arguments:**

    --observer NAME             Observer (e.g., earth, stereo-a, solar-orbiter)
                                Default: model metadata or earth

**Plasma Arguments:**

    --tbase TEMP                Base plasma temperature (K). Default: 1.0e6
    --nbase DENSITY             Base plasma density (cm^-3). Default: 1.0e8
    --q0 FLOAT                  Heating normalization. Default: 0.0217
    --a FLOAT                   Heating exponent a. Default: 0.3
    --b FLOAT                   Heating exponent b. Default: 2.7
    --corona-mode INT           Corona mode bitmask
    --force-isothermal          Enable isothermal mode
    --interpol-b                Enable B-field interpolation
    --analytical-nt             Enable analytical NT mode
    --selective-heating         Enable selective heating
    --shtable-path PATH         Custom 7x7 heating table (.npy, .npz, .txt, .csv)

**Preview & Visualization Arguments:**

    --preview-channel-index INT  Channel index for preview PNG. Default: 0
    --preview-png PATH           Custom preview PNG path
    --log-scale                  Apply log10 scaling to preview
    --omp-threads INT            OpenMP threads for DLL. Default: 8

**Examples:**

    # Basic rendering with automatic defaults
    gxrender-euv --model-path model.h5 --ebtel-path ""

    # With explicit geometry and plasma
    gxrender-euv \
      --model-path model.h5 \
      --ebtel-path /path/to/ebtel.sav \
      --dx 2.0 --dy 2.0 \
      --tbase 1.5e6 --nbase 2.5e8 --q0 0.03 \
      --output-dir /tmp/renders

    # Non-Earth observer with log-scaled preview
    gxrender-euv \
      --model-path model.h5 \
      --observer stereo-a \
      --preview-channel-index 1 \
      --log-scale \
      --output-dir /tmp

    # Multiple channels in preview selection
    gxrender-euv \
      --model-path model.h5 \
      --ebtel-path "" \
      --preview-channel-index 0  # A171
      --preview-png /tmp/preview_a171.png

    # Selective heating with custom table
    gxrender-euv \
      --model-path model.h5 \
      --selective-heating \
      --shtable-path custom_7x7.npy \
      --log-scale

----

Interactive Viewer: gxrender-map-view
--------------------------------------

**Purpose:** View rendered MW/EUV maps with interactive controls.

**Syntax:**

    gxrender-map-view /path/to/maps.h5

**Controls:**

    - Channel slider: Switch between wavelengths (updates colormap automatically)
    - Frequency slider: For MW maps
    - Mouse wheel/trackpad: Zoom
    - Left-click drag: Pan
    - Right-click menu: Color normalization, export
    - 'q' key: Quit viewer

**Features:**

    - Instrument-specific colormaps (AIA, EUVI, TRACE, SXT, EUI)
    - WCS-aware axes with solar coordinates
    - Colorbar with units
    - Interactive legend

----

Environment Variables
---------------------

**GXIMAGECOMPUTING_EBTEL_PATH**

    Path to EBTEL table for DEM/DDM lookups. Used as fallback by CLI.

    export GXIMAGECOMPUTING_EBTEL_PATH="$SSW/packages/gx_simulator/euv/ebtel/ebtel_ss.sav"

**SUNPY_CONFIGDIR**

    SunPy configuration directory (if default is non-writable).

    export SUNPY_CONFIGDIR=/tmp/sunpy_cfg

**MPLCONFIGDIR**

    Matplotlib configuration directory (if default is non-writable).

    export MPLCONFIGDIR=/tmp/mpl_cfg

----

Troubleshooting
---------------

**"Could not find channel A284"**

    Check available channels in the H5 file:
    python -c "import h5py; f = h5py.File('map.h5'); print(list(f['maps']['channel_ids'][:]))"

**Preview looks too dark**

    Try log scaling:
    gxrender-euv-preview map.h5 --log-scale

**"zsh: command not found: gxrender-euv-preview"**

    Reinstall the package:
    pip install -e .

**Matplotlib window won't display**

    Set display environment variable (SSH/headless):
    export MPLBACKEND=Agg
    gxrender-euv-preview map.h5 --no-save

**File size issues with log-scaled PNGs**

    Log-scale PNGs are ~2.5× larger (more color variation).
    This is normal and improves visualization quality.

----

Output Artifacts
----------------

**Rendered H5 Map File:**

    <output-dir>/<model-name>_py_euv_maps.h5

    Contains:
    - /maps/data: [nx, ny, n_channels, 2] flux array (corona, TR)
    - /maps/channel_ids: Channel labels
    - /metadata: WCS header, instrument, observer, dates

**Preview PNG:**

    <output-dir>/<model-name>_py_euv_maps_preview_<CHANNEL>.png

    3-panel figure: Corona | Transition Region | Total
    Uses instrument-specific colormap automatically
    Optional log10 scaling applied if --log-scale specified

----

API Usage (Python)
------------------

**Standalone preview from H5:**

    from gxrender.viz.regenerate_euv_preview import main
    import sys
    sys.argv = ['prog', '/path/to/map.h5', '--channel-id', 'A171', '--log-scale']
    main()

**Workflow-level control:**

    from argparse import Namespace
    from gxrender.workflows.render_euv import run

    args = Namespace(
        model_path="/path/to/model.h5",
        model_format="auto",
        ebtel_path="",
        output_dir="/tmp",
        output_name=None,
        preview_channel_index=1,
        preview_png=None,
        log_scale=True,
        # ... other args
    )
    result = run(args)

**SDK-level (programmatic):**

    from gxrender import render_euv_maps, EUVRenderOptions

    result = render_euv_maps(
        EUVRenderOptions(
            model_path="/path/to/model.h5",
            model_format="sav",
            ebtel_path="",
            observer_name="stereo-a",
            save_outputs=True,
            write_preview=True,  # Generate preview
        )
    )

    preview_path = result.outputs.preview_path  # Path to PNG
    h5_path = result.outputs.h5_path            # Path to H5

----

See Also
--------

- Full documentation: docs/euv_preview.rst
- Examples: examples/python/cli/RenderExampleEUV.py
- API reference: docs/cli_api.rst
