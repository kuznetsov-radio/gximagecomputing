#!/usr/bin/env python

"""
Example: Model Loader with Explicit Metadata Overrides

This example demonstrates how to use explicit parameter overrides (DSun, lonC, b0Sun)
when loading model data, and validates that gximagecomputing correctly applies
these overrides in both H5 and SAV loaders, producing equivalent results.

Purpose:
    Show how users can override model metadata (observer distance, central longitude,
    solar B0 angle) when loading models. This is useful when you want to render a
    model with different viewing geometry than what was stored in the file.

Scientific Context:
    - DSun: Sun-observer distance in cm (affects coordinate scaling)
    - lonC: Central longitude for the view (affects coordinate centering)
    - b0Sun: Solar B0 angle at observation time (affects vertical coordinate mapping)

    These overrides allow renderers to simulate observations from different locations
    or times without needing to regenerate the entire model.

Requirements:
    - pyGXrender-test-data repository available at ../pyGXrender-test-data/raw/
    - Model files with .sav and .clone.h5 variants
"""

import os
import sys
import tempfile
from pathlib import Path

# Configure temporary directories for dependencies
os.environ.setdefault(
    "SUNPY_CONFIGDIR",
    str(Path(tempfile.gettempdir()) / "gximagecomputing_sunpy_config"),
)
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "gximagecomputing_mpl_config"),
)

import numpy as np

from gxrender.io.model import load_model_hdf_with_metadata, load_model_sav_with_metadata
from gxrender.utils.test_data import try_find_model_loader_parity_files


def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def main():
    # Find test data fixtures
    parity_files = try_find_model_loader_parity_files()
    if parity_files is None:
        print("ERROR: Model loader parity test fixtures not found.")
        print("Expected: SAV and .clone.H5 files in pyGXrender-test-data/raw/")
        return False

    sav_path, h5_clone_path = parity_files

    print_section("Model Loader with Explicit Metadata Overrides")

    # Define override values
    override_dsun_cm = 1.4321098765e13
    override_lonc_deg = 23.456789
    override_b0sun_deg = -6.54321

    print(f"\nTest scenario: Load model with custom metadata values:")
    print(f"  DSun override:   {override_dsun_cm:.6e} cm")
    print(f"  lonC override:   {override_lonc_deg:.6f} deg")
    print(f"  b0Sun override:  {override_b0sun_deg:.6f} deg")

    print(f"\nLoading models from:")
    print(f"  SAV file: {sav_path}")
    print(f"  H5 file:  {h5_clone_path}")

    # Load both models with identical overrides
    model_h5, _, metadata_h5 = load_model_hdf_with_metadata(
        str(h5_clone_path),
        DSun=override_dsun_cm,
        lonC=override_lonc_deg,
        b0Sun=override_b0sun_deg,
    )

    model_sav, _, metadata_sav = load_model_sav_with_metadata(
        str(sav_path),
        DSun=override_dsun_cm,
        lonC=override_lonc_deg,
        b0Sun=override_b0sun_deg,
    )

    print_section("1. Verifying Override Application")

    # Check H5 model
    h5_dsun = float(model_h5["DSun"][0])
    h5_lonc = float(model_h5["lonC"][0])
    h5_b0sun = float(model_h5["b0Sun"][0])

    print(f"\nH5 Model Field Values:")
    print(f"  DSun stored in model: {h5_dsun:.6e} cm")
    print(f"  lonC stored in model: {h5_lonc:.6f} deg")
    print(f"  b0Sun stored in model: {h5_b0sun:.6f} deg")

    h5_dsun_match = np.isclose(h5_dsun, override_dsun_cm)
    h5_lonc_match = np.isclose(h5_lonc, override_lonc_deg)
    h5_b0sun_match = np.isclose(h5_b0sun, override_b0sun_deg)

    print(f"\n  {'✓' if h5_dsun_match else '✗'} DSun matches override")
    print(f"  {'✓' if h5_lonc_match else '✗'} lonC matches override")
    print(f"  {'✓' if h5_b0sun_match else '✗'} b0Sun matches override")

    # Check SAV model
    sav_dsun = float(model_sav["DSun"][0])
    sav_lonc = float(model_sav["lonC"][0])
    sav_b0sun = float(model_sav["b0Sun"][0])

    print(f"\nSAV Model Field Values:")
    print(f"  DSun stored in model: {sav_dsun:.6e} cm")
    print(f"  lonC stored in model: {sav_lonc:.6f} deg")
    print(f"  b0Sun stored in model: {sav_b0sun:.6f} deg")

    sav_dsun_match = np.isclose(sav_dsun, override_dsun_cm)
    sav_lonc_match = np.isclose(sav_lonc, override_lonc_deg)
    sav_b0sun_match = np.isclose(sav_b0sun, override_b0sun_deg)

    print(f"\n  {'✓' if sav_dsun_match else '✗'} DSun matches override")
    print(f"  {'✓' if sav_lonc_match else '✗'} lonC matches override")
    print(f"  {'✓' if sav_b0sun_match else '✗'} b0Sun matches override")

    print_section("2. Verifying H5 and SAV Consistency")

    h5_sav_dsun_match = np.isclose(h5_dsun, sav_dsun)
    h5_sav_lonc_match = np.isclose(h5_lonc, sav_lonc)
    h5_sav_b0sun_match = np.isclose(h5_b0sun, sav_b0sun)

    print(f"\nH5 vs SAV comparison (with overrides applied):")
    print(f"  {'✓' if h5_sav_dsun_match else '✗'} DSun matches between formats")
    print(f"  {'✓' if h5_sav_lonc_match else '✗'} lonC matches between formats")
    print(f"  {'✓' if h5_sav_b0sun_match else '✗'} b0Sun matches between formats")

    print_section("3. Verifying Metadata Consistency")

    print(f"\nMetadata values (from structured array parsing):")
    print(f"  H5 metadata DSun:  {metadata_h5['DSun']:.6e} cm")
    print(f"  SAV metadata DSun: {metadata_sav['DSun']:.6e} cm")
    print(f"  {'✓' if metadata_h5['DSun'] == metadata_sav['DSun'] else '✗'} Metadata DSun matches")

    print(f"\n  H5 metadata lonC:  {metadata_h5['lonC']:.6f} deg")
    print(f"  SAV metadata lonC: {metadata_sav['lonC']:.6f} deg")
    print(f"  {'✓' if metadata_h5['lonC'] == metadata_sav['lonC'] else '✗'} Metadata lonC matches")

    print(f"\n  H5 metadata b0Sun:  {metadata_h5['b0Sun']:.6f} deg")
    print(f"  SAV metadata b0Sun: {metadata_sav['b0Sun']:.6f} deg")
    print(f"  {'✓' if metadata_h5['b0Sun'] == metadata_sav['b0Sun'] else '✗'} Metadata b0Sun matches")

    print_section("Summary")

    all_match = (
        h5_dsun_match and h5_lonc_match and h5_b0sun_match and
        sav_dsun_match and sav_lonc_match and sav_b0sun_match and
        h5_sav_dsun_match and h5_sav_lonc_match and h5_sav_b0sun_match and
        metadata_h5['DSun'] == metadata_sav['DSun'] and
        metadata_h5['lonC'] == metadata_sav['lonC'] and
        metadata_h5['b0Sun'] == metadata_sav['b0Sun']
    )

    if all_match:
        print("\n✓ Override parameters correctly applied and consistent across formats!")
        print("\n  This validates that gximagecomputing correctly handles metadata")
        print("  overrides in both H5 and SAV loaders.")
        return True
    else:
        print("\n✗ Inconsistencies found in override application or format comparison.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
