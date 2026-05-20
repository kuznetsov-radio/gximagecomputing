EUV Preview and Colormap Features
==================================

This guide covers the new EUV preview generation and instrument-specific colormap features introduced in pyGXrender v0.0.3.0.

--------

Overview
--------

Three major enhancements improve EUV visualization:

1. **Instrument-Specific Colormaps**: SunPy's 73 predefined EUV colormaps automatically selected based on instrument and channel
2. **Standalone Preview CLI**: New ``gxrender-euv-preview`` command for generating channel-specific previews from rendered H5 files
3. **Log-Scale Rendering**: ``--log-scale`` option reveals the full dynamic range of coronal structures

--------

Instrument-Specific Colormaps
=============================

SunPy Colormap Mapping
----------------------

The viewer (``gxrender-map-view``) and preview generation now automatically select instrument-appropriate colormaps:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Instrument
     - Example Channels
     - SunPy Colormaps Used
   * - AIA
     - 171, 193, 211, 335, 94, 131
     - ``sdoaia171``, ``sdoaia193``, ``sdoaia211``, ``sdoaia335``, ``sdoaia94``, ``sdoaia131``
   * - EUVI (STEREO-A/B)
     - 171, 195, 284, 304
     - ``euvi171``, ``euvi195``, ``euvi284``, ``euvi304``
   * - TRACE
     - 171, 195, 284, 1600
     - ``trace171``, ``trace195``, ``trace284``
   * - SXT (Yohkoh)
     - Al, Al-Mg, Be
     - ``yohkohsxtal``, ``yohkohsxtalmg``, ``yohkohsxtbe``
   * - EUI (Solar Orbiter)
     - FSI-174, HRI-EUV
     - ``solar-orbiter-eui-fsi-174``, ``solar-orbiter-eui-hri-euv``

Fallback
--------

If no matching SunPy colormap is found, the viewer defaults to ``magma`` (a perceptually uniform colormap).

The colormap is automatically updated when:

- Switching channels with the slider in ``gxrender-map-view``
- Generating previews with different channel selections

--------

Preview Generation: Core Workflow
==================================

The ``gxrender-euv`` workflow now supports channel-specific preview PNG generation during rendering:

.. code-block:: bash

   gxrender-euv \
     --model-path model.h5 \
     --ebtel-path "" \
     --preview-channel-index 1 \
     --preview-png /tmp/my_custom_preview.png \
     --log-scale

**Options:**

- ``--preview-channel-index`` (int, default 0): Select which response channel to preview by index
- ``--preview-png`` (Path): Custom output path for preview PNG. Default: ``<model>_py_euv_maps_preview_<channel>.png``
- ``--log-scale`` (flag): Apply log10 scaling for enhanced dynamic range

The preview shows three panels:

1. **Corona** (upper 0.5–3 MK material)
2. **Transition Region** (0.01–0.5 MK material)
3. **Total** (Corona + TR)

Each panel uses the instrument's native SunPy colormap automatically.

--------

Standalone Preview CLI: gxrender-euv-preview
==============================================

For regenerating or creating new previews from existing rendered H5 files without re-running the full workflow:

.. code-block:: bash

   gxrender-euv-preview /path/to/rendered_euv_maps.h5 [OPTIONS]

Channel Selection
-----------------

Select the preview channel in two ways:

**By channel label (recommended):**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A171
   gxrender-euv-preview map.h5 --channel-id 195    # Both "A195" and "195" work

**By index:**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-index 0   # First channel
   gxrender-euv-preview map.h5 --channel-index 2   # Third channel

If no channel is specified, index 0 is used by default.

Output Options
--------------

**Save to custom path:**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A284 --output /tmp/euv_a284.png

**Save to default location (beside H5):**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A284
   # Output: /path/to/map_preview_A284.png

**Display in-memory without saving:**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A171 --no-save --show
   # Matplotlib window opens; press 'q' to close

**Suppress saving, show on screen:**

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A304 --show

Full CLI Help
-------------

.. code-block:: bash

   gxrender-euv-preview --help

   usage: gxrender-euv-preview [-h] [--channel-index CHANNEL_INDEX]
                               [--channel-id CHANNEL_ID] [--output OUTPUT]
                               [--show] [--no-save] [--title TITLE]
                               [--log-scale]
                               map_path

   Regenerate an EUV preview PNG from a rendered H5 map file.

   positional arguments:
     map_path              Rendered EUV map H5 path.

   optional arguments:
     -h, --help            show this help message and exit
     --channel-index CHANNEL_INDEX
                           Channel index to preview (default: 0).
     --channel-id CHANNEL_ID
                           Explicit channel label (for example A171).
     --output OUTPUT       Optional output PNG path. Default:
                           <map-stem>_preview_<channel>.png beside the map.
     --show                Display preview in-memory with matplotlib.
     --no-save             Do not save PNG to disk.
     --title TITLE         Optional preview title override.
     --log-scale           Apply log10 scaling to preview PNG for enhanced
                           dynamic range visualization.

--------

Log-Scale Rendering
====================

The ``--log-scale`` flag applies log₁₀ scaling to dramatically improve visualization of the vast dynamic range in solar EUV data.

Linear vs Log Scale Comparison
------------------------------

**Linear Scale (default):**

- Bright regions dominate the color palette
- Faint structures invisible due to dynamic range compression
- File size: ~95 KB (less color variation)

**Log Scale:**

- Full coronal structure visible
- Transition region details clear
- Faint and bright structures equally visible
- File size: ~250 KB (much more color variation captured)

.. code-block:: bash

   # Linear (default)
   gxrender-euv-preview map.h5 --channel-id A284 --output linear.png

   # Log scale
   gxrender-euv-preview map.h5 --channel-id A284 --output log.png --log-scale

The log scale is applied before creating the SunPy map, ensuring the colormap is applied to the full dynamic range of the transformed data.

Implementation Details
----------------------

Log scaling handles edge cases gracefully:

- NaN and infinite values are masked
- Zero and negative values use the minimum positive value as a floor
- Safe log₁₀ application prevents mathematical errors

This ensures valid visualization even with unusual data distributions.

--------

Usage Examples
==============

Example 1: Render STEREO-A Model with Log-Scaled Preview
---------------------------------------------------------

.. code-block:: bash

   gxrender-euv \
     --model-path stereo_a_model.h5 \
     --ebtel-path /path/to/ebtel_ss.sav \
     --observer stereo-a \
     --preview-channel-index 1 \
     --log-scale \
     --output-dir /tmp/renders

   # Output: /tmp/renders/stereo_a_model.h5_py_euv_maps.h5
   #         /tmp/renders/stereo_a_model.h5_py_euv_maps_preview_A195.png (log scale)

Example 2: Create Multiple Preview Variations from One H5 File
---------------------------------------------------------------

.. code-block:: bash

   # Generate linear previews for all channels
   gxrender-euv-preview map.h5 --channel-id A171 --output preview_a171_linear.png
   gxrender-euv-preview map.h5 --channel-id A195 --output preview_a195_linear.png
   gxrender-euv-preview map.h5 --channel-id A284 --output preview_a284_linear.png
   gxrender-euv-preview map.h5 --channel-id A304 --output preview_a304_linear.png

   # Generate log-scale versions
   gxrender-euv-preview map.h5 --channel-id A171 --output preview_a171_log.png --log-scale
   gxrender-euv-preview map.h5 --channel-id A195 --output preview_a195_log.png --log-scale
   gxrender-euv-preview map.h5 --channel-id A284 --output preview_a284_log.png --log-scale
   gxrender-euv-preview map.h5 --channel-id A304 --output preview_a304_log.png --log-scale

Example 3: Interactive Preview Display
---------------------------------------

.. code-block:: bash

   # View channel A171 without saving to disk
   gxrender-euv-preview map.h5 --channel-id A171 --no-save --show
   # A matplotlib window appears; interact with it (zoom, pan)
   # Press 'q' or close the window to exit

Example 4: Python API Usage
----------------------------

Preview generation can also be called from Python code:

.. code-block:: python

   from pathlib import Path
   from gxrender.workflows.render_euv import _preview_euv
   from astropy.io import fits
   import numpy as np

   # Load your flux data and WCS header
   flux_cor = np.random.rand(512, 512)  # Corona flux
   flux_tr = np.random.rand(512, 512)   # Transition region flux

   # Create a minimal WCS header
   header = fits.Header()
   header['NAXIS'] = 2
   header['NAXIS1'] = 512
   header['NAXIS2'] = 512
   header['DATE-OBS'] = '2020-11-26T19:58:31'

   # Generate preview
   _preview_euv(
       flux_cor=flux_cor,
       flux_tr=flux_tr,
       instrument='EUVIA',
       channel='A171',
       out_png=Path('/tmp/preview.png'),
       title='My EUV Preview',
       base_wcs_header=header,
       observer_text='stereo-a',
       show=False,
       log_scale=True,  # Apply log scaling
   )

--------

Interactive Viewer: gxrender-map-view
======================================

The interactive viewer also benefits from colormap improvements:

.. code-block:: bash

   gxrender-map-view /path/to/rendered_euv_maps.h5

**Features:**

- Channel slider updates both data and colormap dynamically
- Frequency slider (for MW maps)
- Instrument-specific colormaps applied automatically
- Pan, zoom, and color normalization controls

--------

Common Use Cases
================

**Q: How do I compare linear vs log-scale renderings?**

A: Generate both versions and compare side-by-side:

.. code-block:: bash

   gxrender-euv-preview map.h5 --channel-id A284 \
     --output /tmp/linear.png
   gxrender-euv-preview map.h5 --channel-id A284 \
     --output /tmp/log.png --log-scale

**Q: My EUV preview looks too dark. What should I do?**

A: Try log scaling, which reveals faint structures:

.. code-block:: bash

   gxrender-euv-preview map.h5 --log-scale

**Q: Can I change the colormap?**

A: Colormaps are automatically selected based on instrument + channel. To use a different colormap, you can modify the preview generation code or file a feature request.

**Q: How do I know which channel index corresponds to which label?**

A: Check the H5 file's channel_ids dataset:

.. code-block:: python

   import h5py
   with h5py.File('map.h5', 'r') as f:
       print("Available channels:", f['maps']['channel_ids'][:])

**Q: Are the colormap selections accurate to observatories?**

A: Yes, colormaps come from SunPy and are instrument-calibrated. For example:
- AIA colormaps reflect each channel's physical wavelength
- EUVI colormaps match STEREO's ESC instrument
- TRACE colormaps use historical calibration

--------

References
==========

- `SunPy Colormaps Gallery <https://docs.sunpy.org/en/stable/generated/gallery/map/plotting/plot_colormap_selector.html>`_
- `SunPy Map Visualization <https://docs.sunpy.org/en/stable/generated/gallery/map/index.html>`_
- `NumPy Log Scaling <https://numpy.org/doc/stable/reference/generated/numpy.log10.html>`_
