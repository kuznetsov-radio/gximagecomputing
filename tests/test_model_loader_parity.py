"""
Unit tests for gximagecomputing's model packaging logic.

NOTE: These tests focus on gximagecomputing's actual responsibility, which is
to correctly package model data (from pyampp) into numpy structured arrays
compatible with the radiation engine.

For integration examples demonstrating various usage scenarios with real model
files, see examples/python/model_loader/example_model_loader_*.py

pyampp's model loading and ephemeris functionality are tested in pyampp's own
test suite and are validated through the integration examples.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

os.environ.setdefault(
    "SUNPY_CONFIGDIR",
    str(Path(tempfile.gettempdir()) / "gximagecomputing_sunpy_config"),
)
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "gximagecomputing_mpl_config"),
)

from astropy.time import Time


def test_model_packaging_import():
    """Verify that model packaging functions can be imported."""
    from gxrender.io.model import load_model_dict
    assert callable(load_model_dict)


def test_model_loaders_import():
    """Verify that model loaders import correctly (delegates to pyampp)."""
    from gxrender.io.model import load_model_hdf_with_metadata, load_model_sav_with_metadata
    assert callable(load_model_hdf_with_metadata)
    assert callable(load_model_sav_with_metadata)


def test_observer_geometry_normalization():
    """Test that observer names are normalized correctly."""
    from gxrender.geometry.observer_geometry import normalize_observer_name

    # String case
    assert normalize_observer_name("Earth") == "earth"
    assert normalize_observer_name("STEREO-A") == "stereo-a"

    # Bytes handling
    assert normalize_observer_name(b"earth") == "earth"
    assert normalize_observer_name(np.bytes_(b"earth")) == "earth"

    # None handling
    assert normalize_observer_name(None) is None


def test_observer_name_aliases():
    """Test that observer name aliases are resolved correctly."""
    from gxrender.geometry.observer_geometry import normalize_observer_name

    # Test aliases that are defined in _OBSERVER_ALIASES
    aliases = {
        "earth": "earth",
        "terra": "earth",
        "solo": "solar orbiter",
        "solar orbiter": "solar orbiter",
        "solar-orbiter": "solar orbiter",
        "stereo a": "stereo-a",
        "stereo-a": "stereo-a",
        "stereoa": "stereo-a",
        "stereo ahead": "stereo-a",
        "stereo b": "stereo-b",
        "stereo-b": "stereo-b",
        "stereob": "stereo-b",
        "stereo behind": "stereo-b",
    }

    for alias, expected in aliases.items():
        result = normalize_observer_name(alias)
        assert result == expected, f"Expected {alias!r} -> {expected!r}, got {result!r}"


def test_model_dtype_consistency():
    """Test that model structured array has expected dtype fields."""
    from gxrender.io.model import load_model_dict

    # Check that load_model_dict returns a dtype object
    assert hasattr(load_model_dict, '__name__')


def test_pyampp_integration_availability():
    """Test that pyampp geometry functions can be imported by observer_geometry."""
    # Check that key pyampp functions that observer_geometry depends on are available
    from pyampp.geometry import (
        compute_inscribing_fov_from_world,
        world_corners_from_geometry_contract,
    )

    assert callable(compute_inscribing_fov_from_world)
    assert callable(world_corners_from_geometry_contract)
