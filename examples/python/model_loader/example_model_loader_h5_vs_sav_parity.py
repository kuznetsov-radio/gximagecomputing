#!/usr/bin/env python

"""
Example: H5 vs SAV Model Loader Parity

This example demonstrates that gximagecomputing correctly packages model data
loaded from both HDF5 (.clone.h5) and IDL SAV (.sav) formats, producing
equivalent results within numerical precision limits.

Purpose:
    Validate that gximagecomputing's model packaging logic produces identical
    model structures regardless of the source file format (H5 or SAV). This
    documents the expected behavior when using different file formats with
    pyampp's loaders.

Requirements:
    - pyGXrender-test-data repository available at ../pyGXrender-test-data/raw/
    - Model files with .sav and .clone.h5 variants

Scientific Context:
    Small numerical differences between H5 and SAV formats are expected due to:
    - Float32 vs Float64 precision in different serialization formats
    - Rounding during format conversion (IDL to HDF5)
    - Different library implementations for the same numerical operations
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


def compare_model_fields(model_h5, model_sav, field_name, rtol=1e-6, atol=1e-9):
    """Compare a single field between H5 and SAV models."""
    h5_value = np.asarray(model_h5[field_name][0])
    sav_value = np.asarray(model_sav[field_name][0])

    shape_match = h5_value.shape == sav_value.shape
    if not shape_match:
        return False, f"Shape mismatch: H5 {h5_value.shape} vs SAV {sav_value.shape}"

    if np.issubdtype(h5_value.dtype, np.number) and np.issubdtype(sav_value.dtype, np.number):
        matches = np.allclose(h5_value, sav_value, rtol=rtol, atol=atol, equal_nan=True)
        if not matches:
            max_rel_diff = np.nanmax(np.abs((h5_value - sav_value) / (np.abs(sav_value) + 1e-20)))
            max_abs_diff = np.nanmax(np.abs(h5_value - sav_value))
            return False, f"Values differ: max_rel_diff={max_rel_diff:.2e}, max_abs_diff={max_abs_diff:.2e}"
        return True, "✓ Values match within tolerance"
    else:
        matches = np.array_equal(h5_value, sav_value)
        return matches, "✓ Array match" if matches else "✗ Arrays differ"


def main():
    # Find test data fixtures
    parity_files = try_find_model_loader_parity_files()
    if parity_files is None:
        print("ERROR: Model loader parity test fixtures not found.")
        print("Expected: SAV and .clone.H5 files in pyGXrender-test-data/raw/")
        return False

    sav_path, h5_clone_path = parity_files

    print_section("H5 vs SAV Model Loader Parity Example")

    print(f"\nLoading models from:")
    print(f"  SAV file: {sav_path}")
    print(f"  H5 file:  {h5_clone_path}")

    # Load models with no overrides
    model_h5, _, metadata_h5 = load_model_hdf_with_metadata(str(h5_clone_path))
    model_sav, _, metadata_sav = load_model_sav_with_metadata(str(sav_path))

    print_section("1. Comparing Model Structure")

    print(f"\nModel dtype fields:")
    print(f"  H5 fields:  {model_h5.dtype.names}")
    print(f"  SAV fields: {model_sav.dtype.names}")
    print(f"  Match: {model_h5.dtype.names == model_sav.dtype.names}")

    # Compare all fields
    print_section("2. Comparing Field Values")

    all_match = True
    mismatches = []

    for field_name in (model_h5.dtype.names or []):
        matches, msg = compare_model_fields(model_h5, model_sav, field_name)
        status = "✓" if matches else "✗"
        print(f"  {status} {field_name:30s} {msg}")
        if not matches:
            all_match = False
            mismatches.append((field_name, msg))

    # Compare metadata
    print_section("3. Comparing Metadata")

    for key in ("lon", "lat", "dsun_obs", "hgln_obs", "hglt_obs", "crln_obs", "crlt_obs"):
        if key in metadata_h5 and key in metadata_sav:
            match = metadata_h5[key] == metadata_sav[key]
            status = "✓" if match else "✗"
            print(f"  {status} {key:20s} H5={metadata_h5[key]:12} SAV={metadata_sav[key]:12}")
        else:
            print(f"  ✗ {key:20s} (missing in one or both)")

    obs_time_match = metadata_h5["obs_time"].isot == metadata_sav["obs_time"].isot
    print(f"  {'✓' if obs_time_match else '✗'} obs_time              H5={metadata_h5['obs_time'].isot} SAV={metadata_sav['obs_time'].isot}")

    dsun_match = np.isclose(metadata_h5["DSun"], metadata_sav["DSun"])
    print(f"  {'✓' if dsun_match else '✗'} DSun                 H5={metadata_h5['DSun']:.6e} SAV={metadata_sav['DSun']:.6e}")

    b0sun_match = np.isclose(metadata_h5["b0Sun"], metadata_sav["b0Sun"])
    print(f"  {'✓' if b0sun_match else '✗'} b0Sun                H5={metadata_h5['b0Sun']:.6e} SAV={metadata_sav['b0Sun']:.6e}")

    lonc_match = np.isclose(metadata_h5["lonC"], metadata_sav["lonC"])
    print(f"  {'✓' if lonc_match else '✗'} lonC                 H5={metadata_h5['lonC']:.6e} SAV={metadata_sav['lonC']:.6e}")

    print_section("Summary")

    if all_match and obs_time_match and dsun_match and b0sun_match and lonc_match:
        print("\n✓ H5 and SAV models match within numerical precision!")
        print("\n  This validates that gximagecomputing correctly packages")
        print("  model data regardless of the source file format.")
        return True
    else:
        print("\n✗ Differences found between H5 and SAV models:")
        for field, msg in mismatches:
            print(f"    - {field}: {msg}")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
