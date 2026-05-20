#!/usr/bin/env python

"""
Example: Model Loader with Ephemeris Recomputation

This example demonstrates how to use ephemeris recomputation when loading models,
and validates that gximagecomputing correctly applies recomputed observer geometry
in both H5 and SAV loaders.

Purpose:
    Show how users can recompute the observer's ephemeris (position, orientation)
    at the observation time using pyampp's geometry functions. This is useful when
    the stored metadata needs to be refreshed or when simulating observations from
    different observers (Earth, STEREO-A, Solar Orbiter, etc.).

Scientific Context:
    When recompute_observer_ephemeris=True, the loaders call pyampp geometry functions
    to recalculate:
    - Observer position (DSun, observer_lon_deg, observer_lat_deg)
    - Solar orientation angles (b0Sun, lonC)

    This ensures that the viewing geometry is consistent with current ephemeris data,
    rather than relying on potentially outdated values stored in the model file.

Requirements:
    - pyGXrender-test-data repository available at ../pyGXrender-test-data/raw/
    - Model files with .sav and .clone.h5 variants
    - Internet access (for HORIZONS ephemeris queries if cache misses)
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

    print_section("Model Loader with Ephemeris Recomputation (Earth)")

    print(f"\nTest scenario: Load model with recomputed Earth ephemeris")
    print(f"  (Recalculates observer position and orientation from astro data)")

    print(f"\nLoading models from:")
    print(f"  SAV file: {sav_path}")
    print(f"  H5 file:  {h5_clone_path}")

    # Load both models with ephemeris recomputation
    print("\nLoading H5 model with recomputed ephemeris...")
    model_h5, _, metadata_h5 = load_model_hdf_with_metadata(
        str(h5_clone_path),
        recompute_observer_ephemeris=True,
        observer_name="earth",
    )

    print("Loading SAV model with recomputed ephemeris...")
    model_sav, _, metadata_sav = load_model_sav_with_metadata(
        str(sav_path),
        recompute_observer_ephemeris=True,
        observer_name="earth",
    )

    print_section("1. Model Structure Verification")

    print(f"\nModel dtype fields match: {model_h5.dtype.names == model_sav.dtype.names}")
    if model_h5.dtype.names != model_sav.dtype.names:
        print(f"  H5 fields:  {model_h5.dtype.names}")
        print(f"  SAV fields: {model_sav.dtype.names}")
        return False

    print_section("2. Comparing Model Field Values")

    # Compare all numeric fields
    all_match = True
    for field_name in (model_h5.dtype.names or []):
        h5_value = np.asarray(model_h5[field_name][0])
        sav_value = np.asarray(model_sav[field_name][0])

        shape_match = h5_value.shape == sav_value.shape
        if not shape_match:
            print(f"  ✗ {field_name:30s} Shape mismatch: {h5_value.shape} vs {sav_value.shape}")
            all_match = False
            continue

        if np.issubdtype(h5_value.dtype, np.number) and np.issubdtype(sav_value.dtype, np.number):
            # Use 1e-6 relative tolerance (slightly relaxed for recomputed values)
            matches = np.allclose(h5_value, sav_value, rtol=1e-6, atol=1e-6, equal_nan=True)
            status = "✓" if matches else "✗"

            if not matches:
                max_rel_diff = np.nanmax(np.abs((h5_value - sav_value) / (np.abs(sav_value) + 1e-20)))
                max_abs_diff = np.nanmax(np.abs(h5_value - sav_value))
                print(f"  {status} {field_name:30s} max_rel_diff={max_rel_diff:.2e}, max_abs_diff={max_abs_diff:.2e}")
                all_match = False
            else:
                print(f"  {status} {field_name:30s} Values match")
        else:
            matches = np.array_equal(h5_value, sav_value)
            status = "✓" if matches else "✗"
            print(f"  {status} {field_name:30s} {'Array match' if matches else 'Arrays differ'}")
            if not matches:
                all_match = False

    print_section("3. Verifying Geometry Consistency")

    print(f"\nObserver geometry (recomputed):")

    dsun_match = np.isclose(metadata_h5["DSun"], metadata_sav["DSun"])
    print(f"  {'✓' if dsun_match else '✗'} DSun       H5={metadata_h5['DSun']:.6e} cm | SAV={metadata_sav['DSun']:.6e} cm")

    lonc_match = np.isclose(metadata_h5["lonC"], metadata_sav["lonC"])
    print(f"  {'✓' if lonc_match else '✗'} lonC       H5={metadata_h5['lonC']:.6f} deg | SAV={metadata_sav['lonC']:.6f} deg")

    b0sun_match = np.isclose(metadata_h5["b0Sun"], metadata_sav["b0Sun"])
    print(f"  {'✓' if b0sun_match else '✗'} b0Sun      H5={metadata_h5['b0Sun']:.6f} deg | SAV={metadata_sav['b0Sun']:.6f} deg")

    print_section("4. Observation Context")

    print(f"\nObservation time: {metadata_h5['obs_time'].isot}")
    print(f"Observer: Earth (recomputed using pyampp geometry functions)")

    if "observer_lon_deg" in metadata_h5 and "observer_lon_deg" in metadata_sav:
        obs_lon_match = np.isclose(metadata_h5["observer_lon_deg"], metadata_sav["observer_lon_deg"])
        print(f"  Observer longitude: {metadata_h5['observer_lon_deg']:.6f}° (match: {obs_lon_match})")

    if "observer_lat_deg" in metadata_h5 and "observer_lat_deg" in metadata_sav:
        obs_lat_match = np.isclose(metadata_h5["observer_lat_deg"], metadata_sav["observer_lat_deg"])
        print(f"  Observer latitude:  {metadata_h5['observer_lat_deg']:.6f}° (match: {obs_lat_match})")

    if "observer_dsun_cm" in metadata_h5 and "observer_dsun_cm" in metadata_sav:
        obs_dsun_match = np.isclose(metadata_h5["observer_dsun_cm"], metadata_sav["observer_dsun_cm"])
        print(f"  Observer distance:  {metadata_h5['observer_dsun_cm']:.6e} cm (match: {obs_dsun_match})")

    print_section("Summary")

    if all_match and dsun_match and lonc_match and b0sun_match:
        print("\n✓ Ephemeris recomputation produces consistent results across formats!")
        print("\n  This validates that gximagecomputing correctly handles ephemeris")
        print("  recomputation in both H5 and SAV loaders, ensuring that the")
        print("  recomputed observer geometry is consistently applied.")
        return True
    else:
        print("\n✗ Inconsistencies found in ephemeris recomputation or field comparison.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
