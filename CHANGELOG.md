# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- EUV rendering accepts `--parallel`, `--exact`, and `--projection-threads`,
  also on `EUVRenderOptions`. Defaults stay off (`parallel=False`,
  `exact=False`, `projection_threads=0`). The flags are passed to the renderer
  and stored on the result and in the HDF5 metadata.

### Fixed

- The implicit Python-native AIA response now uses the SolarSoft default
  correction state, `evenorm_chiantifix`, rather than `evenorm` alone. The
  normal EUV workflow and ComputeEUV input dumper pin this state explicitly.
- Default AIA response tests now require all seven export channels, with both
  EVE normalization and the CHIANTI correction enabled.

## [0.0.3.0] - 2026-05-20

### Added

#### EUV Preview and Colormap Features

- **Instrument-Specific Colormaps**: Automatic SunPy colormap selection based on instrument + channel
  - AIA: 6 wavelengths (94–335 Å)
  - STEREO EUVI: 4 wavelengths (171–304 Å)
  - TRACE: 3 wavelengths (171–284 Å)
  - Yohkoh SXT: 3 wavelengths (Al, Al-Mg, Be)
  - Solar Orbiter EUI: 2 modes (FSI-174, HRI-EUV)
  - Fallback to `magma` colormap if no match found

- **New CLI: `gxrender-euv-preview`**
  - Standalone preview regeneration from rendered H5 files
  - Channel selection by label (`--channel-id A171`) or index (`--channel-index 0`)
  - Custom output paths (`--output`)
  - Interactive display (`--show`)
  - In-memory preview without file writing (`--no-save`)

- **Log-Scale Rendering (`--log-scale`)**
  - Applies log₁₀ transformation to enhance dynamic range
  - Reveals faint coronal structures invisible in linear scaling
  - Available in both `gxrender-euv` workflow and `gxrender-euv-preview` CLI
  - Safe handling of NaN, infinity, and non-positive values

- **Interactive Viewer Improvements**
  - `gxrender-map-view` now updates colormaps when switching channels
  - Dynamic colormap selection maintains instrument-specific calibration across all channels

- **Preview Generation in Workflows**
  - `gxrender-euv` workflow now supports:
    - `--preview-channel-index`: Select which channel to preview by index
    - `--preview-png`: Custom preview output path
    - `--log-scale`: Apply log scaling to workflow-generated preview
  - Default preview filename includes channel label: `<model>_py_euv_maps_preview_<CHANNEL>.png`

### Enhanced

- **Native EUV Renderer (`source/EUVmain.cpp`)**
  - Updated the C++ transition-region DEM scaling path used by the native `RenderGRFF` library.
  - Guarded TR angular scaling against non-finite values from zero magnetic-field magnitude, zero LOS denominator, and out-of-range `acos` arguments caused by floating-point drift.
  - Sanitized TR DEM samples and TR response integration terms so non-finite intermediate values produce zero TR contribution instead of propagating NaN/Inf into EUV map outputs.
  - Recompiled and added the matching macOS arm64 runtime binary: `binaries/RenderGRFF_arm64.so`.
  - Reviewer note: this is an intentional native-source change, not just Python wrapper work. The touched code block is explicitly marked with comments in `source/EUVmain.cpp` near the TR angular scaling and TR response integration logic.

- **Policy System**: Extended to support non-Earth observers with time-dependent response providers
  - `script_policy.py`: Checks provider availability before rejecting non-AIA instruments
  - `euv_response_policy.py`: Consults provider registry for implicit time-dependent responses
  - `response_providers.py`: Synthesizes providers from pyEUVTools builders

- **EUV Response Providers**: Integrated pyEUVTools response builders
  - Auto-discovery of time-dependent providers (EUVI, EUI, TRACE, SXT)
  - Implicit provider routing in resolution pipeline
  - Priority-based provider selection

- **Documentation**
  - New comprehensive guide: `docs/euv_preview.rst`
  - Quick reference CLI guide: `docs/euv_cli_quickref.md`
  - Updated main README with feature overview
  - Sphinx documentation index updated with new preview section

### Technical Details

- **Colormap Mapping**: `_euv_sunpy_cmap_name(instrument, channel)` function
  - Case-insensitive instrument matching
  - Channel normalization (strips "A" prefix if present)
  - Returns SunPy colormap name or None (falls back to magma)

- **Log Scaling**: `_apply_log_scaling(flux)` function
  - Masks invalid values (NaN, infinity)
  - Uses minimum positive value as floor for non-positive data
  - Safe log₁₀ application prevents mathematical errors

- **H5 Structure**: Extended metadata in rendered EUV maps
  - `maps/channel_ids`: Channel labels for preview selection
  - `maps/component_ids`: Component identifiers (CORONA, TR)
  - `metadata/instrument`: Instrument name for colormap selection

### Testing

- Rebuilt the native library after the `source/EUVmain.cpp` update and included the macOS arm64 binary artifact for review with the source change.
- Added test cases for policy parity after provider registry changes
- Validated provider discovery for non-Earth observers
- Tested colormap selection across all supported instruments
- Verified log scaling with edge cases (NaN, negative values, zero)
- Confirmed H5 file generation with proper channel metadata

### Bug Fixes

- **Non-Earth EUV Rendering**: Fixed stereo-a/stereo-b/solo-fsi/solo-hri/trace/sxt rendering
  - Issue: "requires explicit RESPONSE_SAV" gate hardcoded non-Earth rejection
  - Solution: Check provider registry before rejecting; accept if provider found

- **Provider Discovery**: Corrected API assumptions about pyEUVTools
  - Issue: Looking for non-existent `get_time_dependent_response_providers()` hook
  - Solution: Synthesize providers from actual builder functions (build_euvi_*, build_eui_*, etc.)

- **Preview Channel Selection**: Channel index now properly clipped to valid range
  - Prevents out-of-bounds errors when requesting invalid channel indices

### Compatibility

- **Backward Compatible**: All new features are optional
  - Existing `gxrender-euv` calls work unchanged
  - Default preview behavior (index 0) maintained if `--preview-channel-index` omitted
  - Linear scaling remains default if `--log-scale` not specified

- **SDK Compatibility**: New parameters optional in programmatic usage
  - `_preview_euv(..., log_scale=False)` has sensible default
  - CLI argument parsing handles missing args gracefully

### Performance

- **File Size Impact**: Log-scaled PNGs ~2.5× larger due to enhanced color variation
  - Linear preview: ~95 KB
  - Log-scaled preview: ~240 KB
  - Larger size reflects richer information content

- **Computation**: Log scaling adds negligible overhead (~10 ms per preview)
  - Single numpy operation per flux array
  - Performed before SunPy map creation

### Documentation URLs

- Main guide: [EUV Preview and Colormap Features](docs/euv_preview.rst)
- Quick reference: [CLI Quick Reference](docs/euv_cli_quickref.md)
- Examples: `examples/python/cli/RenderExampleEUV.py`

### Migration Guide

**For existing users:**

No action required. Existing scripts continue to work unchanged.

**To use new features:**

1. Update package: `pip install -e .` or `pip install pyGXrender --upgrade`
2. Try log scaling: Add `--log-scale` to any `gxrender-euv` or `gxrender-euv-preview` command
3. Select channels: Use `--channel-id A171` in preview generation
4. Interactive viewer: Launch `gxrender-map-view map.h5` to see automatic colormaps

### Known Limitations

- Colormap selection is automatic; custom colormaps require code modification
- Log scaling uses base-10; other logarithms would need API extension
- Preview generation requires H5 file with `maps/channel_ids` dataset (all new renders have this)

### Contributors

- Auto-discovery and provider bridging: EUV response policy system
- Log scaling implementation: Safe numerical handling for edge cases
- Colormap integration: SunPy's extensive instrument calibration database

---

## [0.0.2.3 and earlier] - Previous versions

See git history for older releases.
