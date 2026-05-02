from __future__ import annotations

import os
import tempfile
from pathlib import Path

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

from gxrender.io.model import load_model_hdf_with_metadata, load_model_sav_with_metadata
from gxrender.utils.test_data import test_data_setup_hint as _test_data_setup_hint, try_find_model_loader_parity_files


MODEL_LOADER_PARITY_FILES = try_find_model_loader_parity_files()
if MODEL_LOADER_PARITY_FILES is None:
    pytest.skip(_test_data_setup_hint("loader parity fixtures"), allow_module_level=True)
SAV_PATH, H5_CLONE_PATH = MODEL_LOADER_PARITY_FILES

OVERRIDE_DSUN_CM = 1.4321098765e13
OVERRIDE_LONC_DEG = 23.456789
OVERRIDE_B0SUN_DEG = -6.54321


def _assert_models_equal(model_h5: np.ndarray, model_sav: np.ndarray) -> None:
    assert model_h5.dtype.names == model_sav.dtype.names

    for name in model_h5.dtype.names or ():
        h5_value = np.asarray(model_h5[name][0])
        sav_value = np.asarray(model_sav[name][0])

        assert h5_value.shape == sav_value.shape, name
        if np.issubdtype(h5_value.dtype, np.number) and np.issubdtype(sav_value.dtype, np.number):
            assert np.allclose(h5_value, sav_value, rtol=0.0, atol=0.0, equal_nan=True), name
        else:
            assert np.array_equal(h5_value, sav_value), name


def test_clone_h5_loader_matches_sav_loader_exactly() -> None:
    model_h5, _model_dt_h5, metadata_h5 = load_model_hdf_with_metadata(str(H5_CLONE_PATH))
    model_sav, _model_dt_sav, metadata_sav = load_model_sav_with_metadata(str(SAV_PATH))

    _assert_models_equal(model_h5, model_sav)

    for key in ("lon", "lat", "dsun_obs", "hgln_obs", "hglt_obs", "crln_obs", "crlt_obs"):
        assert key in metadata_h5
        assert metadata_h5[key] == metadata_sav[key]

    assert metadata_h5["obs_time"].isot == metadata_sav["obs_time"].isot
    assert np.isclose(metadata_h5["DSun"], metadata_sav["DSun"])
    assert np.isclose(metadata_h5["b0Sun"], metadata_sav["b0Sun"])
    assert np.isclose(metadata_h5["lonC"], metadata_sav["lonC"])


def test_clone_h5_loader_matches_sav_loader_with_explicit_overrides() -> None:
    model_h5, _model_dt_h5, metadata_h5 = load_model_hdf_with_metadata(
        str(H5_CLONE_PATH),
        DSun=OVERRIDE_DSUN_CM,
        lonC=OVERRIDE_LONC_DEG,
        b0Sun=OVERRIDE_B0SUN_DEG,
    )
    model_sav, _model_dt_sav, metadata_sav = load_model_sav_with_metadata(
        str(SAV_PATH),
        DSun=OVERRIDE_DSUN_CM,
        lonC=OVERRIDE_LONC_DEG,
        b0Sun=OVERRIDE_B0SUN_DEG,
    )

    _assert_models_equal(model_h5, model_sav)

    assert np.isclose(float(model_h5["DSun"][0]), OVERRIDE_DSUN_CM)
    assert np.isclose(float(model_h5["lonC"][0]), OVERRIDE_LONC_DEG)
    assert np.isclose(float(model_h5["b0Sun"][0]), OVERRIDE_B0SUN_DEG)
    assert np.isclose(float(model_sav["DSun"][0]), OVERRIDE_DSUN_CM)
    assert np.isclose(float(model_sav["lonC"][0]), OVERRIDE_LONC_DEG)
    assert np.isclose(float(model_sav["b0Sun"][0]), OVERRIDE_B0SUN_DEG)

    assert metadata_h5["DSun"] == metadata_sav["DSun"] == OVERRIDE_DSUN_CM
    assert metadata_h5["lonC"] == metadata_sav["lonC"] == OVERRIDE_LONC_DEG
    assert metadata_h5["b0Sun"] == metadata_sav["b0Sun"] == OVERRIDE_B0SUN_DEG


def test_clone_h5_loader_matches_sav_loader_with_recomputed_earth_ephemeris() -> None:
    model_h5, _model_dt_h5, metadata_h5 = load_model_hdf_with_metadata(
        str(H5_CLONE_PATH),
        recompute_observer_ephemeris=True,
        observer_name="earth",
    )
    model_sav, _model_dt_sav, metadata_sav = load_model_sav_with_metadata(
        str(SAV_PATH),
        recompute_observer_ephemeris=True,
        observer_name="earth",
    )

    assert model_h5.dtype.names == model_sav.dtype.names
    for name in model_h5.dtype.names or ():
        h5_value = np.asarray(model_h5[name][0])
        sav_value = np.asarray(model_sav[name][0])
        assert h5_value.shape == sav_value.shape, name
        if np.issubdtype(h5_value.dtype, np.number) and np.issubdtype(sav_value.dtype, np.number):
            assert np.allclose(h5_value, sav_value, rtol=1e-6, atol=1e-6, equal_nan=True), name
        else:
            assert np.array_equal(h5_value, sav_value), name
    assert np.isclose(metadata_h5["DSun"], metadata_sav["DSun"])
    assert np.isclose(metadata_h5["lonC"], metadata_sav["lonC"])
    assert np.isclose(metadata_h5["b0Sun"], metadata_sav["b0Sun"])
