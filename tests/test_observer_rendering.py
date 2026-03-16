from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import h5py
import numpy as np
import pytest
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import gxrender.io.model as model_io
from gxrender.geometry import observer_geometry
from gxrender.geometry.observer_geometry import (
    build_ephemeris_from_pb0r,
    build_pb0r_from_ephemeris,
    compute_sunpy_wcs_header,
    resolve_observer_geometry,
)
from gxrender.io.model import load_model_hdf_with_observer
from gxrender.io.maps_h5 import save_h5_maps
from gxrender.sdk import CoronalPlasmaParameters, MWRenderOptions, render_mw_maps
from gxrender.utils.test_data import test_data_setup_hint as _test_data_setup_hint, try_find_model_file
from gxrender.utils.render_map_view import _read_render_h5
from gxrender.workflows import render_mw as render_mw_workflow
from gxrender.workflows._render_common import (
    load_model_and_fov,
    resolve_mw_frequencies,
    resolve_plasma_parameters,
)

TEST_CHR_H5 = try_find_model_file("test.chr.h5")


def _model_stub(obs_unix: float = 0.0) -> np.ndarray:
    model = np.zeros(
        1,
        dtype=[
            ("obstime", np.float64),
            ("DSun", np.float64),
            ("lonC", np.float64),
            ("b0Sun", np.float64),
        ],
    )
    model["obstime"] = obs_unix
    model["DSun"] = 1.495978707e13
    model["lonC"] = 12.5
    model["b0Sun"] = -4.25
    return model


def _args(**overrides) -> argparse.Namespace:
    base = {
        "observer": None,
        "lonc_deg": None,
        "b0sun_deg": None,
        "dsun_cm": None,
        "auto_fov": False,
        "use_saved_fov": False,
    }
    base.update(overrides)
    return argparse.Namespace(**base)


def _require_test_chr_h5() -> Path:
    if TEST_CHR_H5 is None:
        pytest.skip(_test_data_setup_hint("test.chr.h5"))
    return TEST_CHR_H5


def test_observer_resolution_priority_prefers_cli_observer() -> None:
    model = _model_stub()
    metadata = {
        "observer": "mars",
        "lon": 200.0,
        "lat": 10.0,
        "crln_obs": 1.0,
        "crlt_obs": 2.0,
        "dsun_obs": 3.0e11,
    }
    resolved = resolve_observer_geometry(model, _args(observer="earth"), metadata)
    assert resolved.observer_name == "earth"
    assert resolved.observer_source == "cli_observer"
    assert np.isfinite(resolved.dsun_cm)
    assert resolved.render_b0_deg == resolved.b0_deg
    assert resolved.render_dsun_cm == resolved.dsun_cm


def test_observer_resolution_uses_non_earth_body() -> None:
    model = _model_stub()
    resolved = resolve_observer_geometry(model, _args(observer="mars"), {})
    assert resolved.observer_name == "mars"
    assert resolved.observer_source == "cli_observer"
    assert np.isfinite(resolved.l0_deg)
    assert np.isfinite(resolved.b0_deg)
    assert resolved.dsun_cm > 0


def test_spacecraft_observer_falls_back_to_horizons() -> None:
    model = _model_stub()
    original_body = observer_geometry.get_body_heliographic_stonyhurst
    original_horizons = observer_geometry.get_horizons_coord
    try:
        def _raise_body(*_args, **_kwargs):
            raise RuntimeError("no built-in ephemeris")

        def _fake_horizons(body, time, id_type=None, *, include_velocity=False):
            assert body == "Solar Orbiter"
            return SkyCoord(
                lon=15.0 * u.deg,
                lat=-3.0 * u.deg,
                radius=0.8 * u.AU,
                frame="heliographic_stonyhurst",
                obstime=time,
            )

        observer_geometry.get_body_heliographic_stonyhurst = _raise_body
        observer_geometry.get_horizons_coord = _fake_horizons
        resolved = resolve_observer_geometry(model, _args(observer="solo"), {"lon": 200.0, "lat": 10.0})
    finally:
        observer_geometry.get_body_heliographic_stonyhurst = original_body
        observer_geometry.get_horizons_coord = original_horizons

    assert resolved.observer_name == "solar orbiter"
    assert resolved.observer_source == "cli_observer"
    assert resolved.dsun_cm > 0


def test_default_path_uses_metadata_observer_for_render_triad() -> None:
    model = _model_stub()
    metadata = {
        "lon": 27.8,
        "crln_obs": 44.9,
        "crlt_obs": 1.4,
        "dsun_obs": 1.476e11,
    }
    resolved = resolve_observer_geometry(model, _args(), metadata)
    assert resolved.observer_name == "custom"
    assert resolved.observer_source == "model_metadata_carrington"
    assert np.isfinite(resolved.render_lonc_deg)
    assert resolved.render_lonc_deg != float(model["lonC"][0])
    assert np.isclose(resolved.render_b0_deg, resolved.b0_deg)
    assert np.isclose(resolved.render_dsun_cm, resolved.dsun_cm)
    assert np.isfinite(resolved.l0_deg)
    assert np.isfinite(resolved.b0_deg)
    assert not np.isclose(resolved.l0_deg, 44.9)


def test_saved_custom_observer_metadata_is_restored_from_pb0r() -> None:
    model = _model_stub()
    metadata = {
        "observer_name": "custom",
        "observer_b0_deg": 5.0,
        "observer_l0_deg": 30.0,
        "observer_rsun_arcsec": 1000.0,
        "observer_rsun_cm": 6.957e10,
    }
    resolved = resolve_observer_geometry(model, _args(), metadata)
    expected_dsun = 6.957e10 / np.sin((1000.0 * u.arcsec).to_value(u.rad))
    assert resolved.observer_name == "custom"
    assert resolved.observer_source == "saved_observer_metadata"
    assert np.isclose(resolved.b0_deg, 5.0)
    assert np.isclose(resolved.l0_deg, 30.0)
    assert np.isclose(resolved.dsun_cm, expected_dsun)
    assert np.isclose(resolved.rsun_arcsec, 1000.0)


def test_pb0r_ephemeris_roundtrip() -> None:
    ephemeris = build_ephemeris_from_pb0r(
        b0_deg=5.0,
        l0_deg=30.0,
        rsun_arcsec=1000.0,
        obs_date="2025-11-26T15:34:31.400",
        rsun_cm=6.957e10,
    )
    assert ephemeris is not None
    pb0r = build_pb0r_from_ephemeris(ephemeris)
    assert pb0r is not None
    assert np.isclose(pb0r["b0_deg"], 5.0)
    assert np.isclose(pb0r["l0_deg"], 30.0)
    assert np.isclose(pb0r["rsun_arcsec"], 1000.0)


def test_h5_stores_sunpy_compatible_wcs_and_viewer_reads_it(tmp_path) -> None:
    model = _model_stub()
    geometry = resolve_observer_geometry(model, _args(observer="earth"), {})
    wcs_header = compute_sunpy_wcs_header(
        nx=6,
        ny=4,
        xc_arcsec=10.0,
        yc_arcsec=-20.0,
        dx_arcsec=2.4,
        dy_arcsec=2.4,
        obs_time="2024-05-12T16:00:00",
        observer_geometry=geometry,
        bunit="K",
    )
    result = {
        "TI": np.ones((4, 6, 2), dtype=np.float32),
        "TV": np.zeros((4, 6, 2), dtype=np.float32),
    }
    out_h5 = tmp_path / "mw.h5"
    save_h5_maps(
        result=result,
        freqlist=[5.8, 6.0],
        out_h5=out_h5,
        model_path=tmp_path / "model.h5",
        model_format="h5",
        xc=10.0,
        yc=-20.0,
        dx=2.4,
        dy=2.4,
        obs_time_iso="2024-05-12T16:00:00",
        wcs_header=wcs_header,
        observer_name=geometry.observer_name,
        observer_source=geometry.observer_source,
        observer_warnings=geometry.warnings,
        l0_deg=geometry.l0_deg,
        b0_deg=geometry.b0_deg,
        dsun_cm=geometry.dsun_cm,
        rsun_cm=geometry.rsun_cm,
        rsun_arcsec=geometry.rsun_arcsec,
        tbase=2.5e6,
        nbase=3.5e8,
        q0=0.031,
        a=0.45,
        b=2.2,
        corona_mode=5,
        shtable=np.full((7, 7), 0.75, dtype=np.float64),
    )

    with h5py.File(out_h5, "r") as h5:
        assert "wcs_header" in h5["metadata"]
        assert "observer_rsun_arcsec" in h5["metadata"]
        assert np.isclose(h5["metadata"]["tbase_k"][()], 2.5e6)
        assert np.isclose(h5["metadata"]["nbase_cm3"][()], 3.5e8)
        assert np.isclose(h5["metadata"]["q0"][()], 0.031)
        assert np.isclose(h5["metadata"]["a"][()], 0.45)
        assert np.isclose(h5["metadata"]["b"][()], 2.2)
        assert int(h5["metadata"]["corona_mode"][()]) == 5
        assert np.allclose(h5["metadata"]["shtable"][()], 0.75)
        header = fits.Header.fromstring(h5["metadata"]["wcs_header"][()].decode("utf-8"), sep="\n")
        assert header["CTYPE1"] == "HPLN-TAN"
        assert header["CTYPE2"] == "HPLT-TAN"
        assert header["OBSERVER"] == "Earth"
        assert np.isclose(header["B0"], geometry.b0_deg)
        assert np.isclose(header["L0"], geometry.l0_deg)
        assert np.isclose(header["RSUN_ARC"], geometry.rsun_arcsec)
        assert np.isclose(header["SOLAR_B0"], geometry.b0_deg)
        assert np.isclose(header["SOLAR_L0"], geometry.l0_deg)
        assert np.isclose(header["HGLT_OBS"], geometry.b0_deg)
        assert np.isclose(header["HGLN_OBS"], geometry.l0_deg)
        assert np.isclose(header["RSUN_OBS"], geometry.rsun_arcsec)

    data = np.ones((4, 6), dtype=np.float32)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = sunpy.map.Map(data, header)
    assert m.coordinate_frame is not None
    assert len(caught) == 0

    viewer_data = _read_render_h5(out_h5)
    assert viewer_data["index_header"]
    assert viewer_data["observer_name"] == "earth"


def test_saved_observer_fov_is_used_when_present() -> None:
    model_path = _require_test_chr_h5()
    (
        model,
        _model_dt,
        _metadata,
        geometry,
        center_source,
        xc_auto,
        yc_auto,
        model_w_arcsec,
        model_h_arcsec,
        applied,
    ) = load_model_and_fov(model_path, "h5", _args())

    assert geometry.observer_name == "stereo-a"
    assert geometry.observer_source == "saved_observer_metadata"
    assert center_source == "saved_observer_fov"
    assert np.isclose(applied["lonC_deg"], geometry.render_lonc_deg)
    assert np.isclose(applied["b0Sun_deg"], geometry.b0_deg)
    assert np.isclose(model["lonC"][0], geometry.render_lonc_deg)
    assert np.isclose(model["b0Sun"][0], geometry.b0_deg)
    assert np.isclose(xc_auto, -903.0)
    assert np.isclose(yc_auto, -171.0)
    assert np.isclose(model_w_arcsec, 150.0)
    assert np.isclose(model_h_arcsec, 150.0)
    assert np.isclose(applied["DSun_cm"], geometry.render_dsun_cm)


def test_h5_loader_obs_time_ignores_metadata_execute_time() -> None:
    model_path = _require_test_chr_h5()
    model, _model_dt, metadata, _observer = load_model_hdf_with_observer(str(model_path))

    assert metadata["obs_time"].isot == "2025-11-26T15:34:31.400"
    assert metadata["observer_obs_date"] == "2025-11-26T15:34:31.400"
    assert metadata["observer_pb0r_obs_date"] == "2025-11-26T15:34:31.400"
    assert metadata["execute"].startswith("gx-fov2box --time 2025-11-26T15:47:52")
    assert model["obstime"][0] != 1480175272.0


def test_base_index_header_is_read_as_fits_cards(tmp_path) -> None:
    hdr = fits.Header()
    hdr["CRVAL1"] = 27.812563081166793
    hdr["CRVAL2"] = -12.247129818437514
    hdr["DATE-OBS"] = "2025-11-26T15:34:31.400"
    hdr["DSUN_OBS"] = 147639256672.79785
    hdr["OBSERVER"] = "STEREO-A"
    hdr["HGLN_OBS"] = 49.91244944046865
    hdr["HGLT_OBS"] = -4.657776399560966
    hdr["CRLN_OBS"] = 94.66120123875693
    hdr["CRLT_OBS"] = -4.657776399560966

    path = tmp_path / "base_index_only.h5"
    with h5py.File(path, "w") as h5:
        base = h5.create_group("base")
        base.create_dataset("index", data=np.bytes_(hdr.tostring(sep="\n", endcard=False, padding=False)))

    extracted: dict[str, object] = {}
    with h5py.File(path, "r") as h5:
        model_io._fill_header_from_base_index(h5, extracted)

    assert np.isclose(float(extracted["lon"]), hdr["CRVAL1"])
    assert np.isclose(float(extracted["lat"]), hdr["CRVAL2"])
    assert extracted["obs_time"] == hdr["DATE-OBS"]
    assert np.isclose(float(extracted["dsun_obs"]), hdr["DSUN_OBS"])
    assert extracted["observer"] == hdr["OBSERVER"]
    assert np.isclose(float(extracted["hgln_obs"]), hdr["HGLN_OBS"])
    assert np.isclose(float(extracted["hglt_obs"]), hdr["HGLT_OBS"])
    assert np.isclose(float(extracted["crln_obs"]), hdr["CRLN_OBS"])
    assert np.isclose(float(extracted["crlt_obs"]), hdr["CRLT_OBS"])


def test_auto_fov_overrides_saved_observer_fov() -> None:
    model_path = _require_test_chr_h5()
    (
        _model,
        _model_dt,
        _metadata,
        geometry,
        center_source,
        xc_auto,
        yc_auto,
        model_w_arcsec,
        model_h_arcsec,
        _applied,
    ) = load_model_and_fov(model_path, "h5", _args(auto_fov=True))

    assert geometry.observer_name == "stereo-a"
    assert center_source == "inscribing_fov"
    assert np.isclose(xc_auto, -983.70, atol=3.0)
    assert np.isclose(yc_auto, -197.37, atol=3.0)
    assert np.isclose(model_w_arcsec, 335.31, atol=4.0)
    assert np.isclose(model_h_arcsec, 335.31, atol=4.0)


def test_computed_fov_uses_execute_geometry_anchor() -> None:
    model_path = _require_test_chr_h5()
    model, _model_dt, metadata, observer_metadata = load_model_hdf_with_observer(str(model_path))
    geometry = resolve_observer_geometry(model, _args(), metadata, observer_metadata)
    fov = observer_geometry.compute_inscribing_fov(
        model,
        geometry,
        model_metadata=metadata,
        observer_metadata=observer_metadata,
    )

    assert np.isclose(fov["xc_arcsec"], -983.70, atol=3.0)
    assert np.isclose(fov["yc_arcsec"], -197.37, atol=3.0)
    assert np.isclose(fov["xsize_arcsec"], 335.31, atol=4.0)
    assert np.isclose(fov["ysize_arcsec"], 335.31, atol=4.0)


def test_resolve_plasma_and_frequency_overrides() -> None:
    custom_shtable = np.arange(49, dtype=np.float64).reshape(7, 7)
    args = argparse.Namespace(
        frequencies_ghz=[6.2, 7.4, 9.8],
        tbase=2.2e6,
        nbase=4.1e8,
        q0=0.019,
        a=0.5,
        b=2.1,
        corona_mode=2,
        selective_heating=False,
        force_isothermal=True,
        interpol_b=False,
        analytical_nt=True,
        shtable=custom_shtable,
        shtable_path=None,
    )

    assert resolve_mw_frequencies(args) == [6.2, 7.4, 9.8]
    plasma = resolve_plasma_parameters(args)
    assert np.isclose(plasma.tbase, 2.2e6)
    assert np.isclose(plasma.nbase, 4.1e8)
    assert np.isclose(plasma.q0, 0.019)
    assert np.isclose(plasma.a, 0.5)
    assert np.isclose(plasma.b, 2.1)
    assert plasma.mode == 7
    assert np.array_equal(plasma.shtable, custom_shtable)


def test_shared_workflow_requires_explicit_science_inputs() -> None:
    with pytest.raises(ValueError, match="Microwave frequency list must be explicit"):
        resolve_mw_frequencies(argparse.Namespace(frequencies_ghz=None))

    with pytest.raises(ValueError, match="Coronal plasma inputs must be explicit"):
        resolve_plasma_parameters(
            argparse.Namespace(
                tbase=None,
                nbase=4.1e8,
                q0=0.019,
                a=0.5,
                b=2.1,
                corona_mode=None,
                selective_heating=False,
                force_isothermal=False,
                interpol_b=False,
                analytical_nt=False,
                shtable=None,
                shtable_path=None,
            )
        )


def test_shared_workflow_keeps_shtable_undefined_when_not_supplied() -> None:
    plasma = resolve_plasma_parameters(
        argparse.Namespace(
            tbase=1.5e6,
            nbase=2.5e8,
            q0=0.02,
            a=0.4,
            b=2.2,
            corona_mode=None,
            selective_heating=False,
            force_isothermal=False,
            interpol_b=False,
            analytical_nt=False,
            shtable=None,
            shtable_path=None,
        )
    )

    assert plasma.mode == 0
    assert plasma.selective_heating is False
    assert plasma.shtable is None


def test_selective_heating_without_explicit_shtable_uses_default_table() -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        plasma = resolve_plasma_parameters(
            argparse.Namespace(
                tbase=1.5e6,
                nbase=2.5e8,
                q0=0.02,
                a=0.4,
                b=2.2,
                corona_mode=None,
                selective_heating=True,
                force_isothermal=False,
                interpol_b=False,
                analytical_nt=False,
                shtable=None,
                shtable_path=None,
            )
        )

    assert plasma.selective_heating is True
    assert plasma.shtable is not None
    assert np.asarray(plasma.shtable).shape == (7, 7)
    assert any("Selective heating was enabled without an explicit SHtable" in str(w.message) for w in caught)


def test_sdk_mw_options_forward_frequency_and_plasma_controls(monkeypatch, tmp_path) -> None:
    captured = {}
    custom_shtable = np.arange(49, dtype=np.float64).reshape(7, 7)

    def _fake_run(ns, verbose=False):
        captured["ns"] = ns
        captured["verbose"] = verbose
        return {
            "library_path": "libRenderGRFF.so",
            "model_path": str(ns.model_path),
            "model_format": str(ns.model_format),
            "ebtel_path": str(ns.ebtel_path or ""),
            "observer_overrides_applied": {},
            "center_source": "saved_observer_fov",
            "geometry": {
                "xc_arcsec": 10.0,
                "yc_arcsec": -20.0,
                "dx_arcsec": 2.0,
                "dy_arcsec": 2.0,
                "nx": 8,
                "ny": 6,
                "fov_x_arcsec": 16.0,
                "fov_y_arcsec": 12.0,
            },
            "obs_time_iso": "2025-11-26T15:34:31.400",
            "freqlist_ghz": [float(v) for v in ns.frequencies_ghz],
            "plasma": {
                "tbase_k": float(ns.tbase),
                "nbase_cm3": float(ns.nbase),
                "q0": float(ns.q0),
                "a": float(ns.a),
                "b": float(ns.b),
                "mode": int(ns.corona_mode),
                "shtable": np.asarray(ns.shtable, dtype=np.float64).tolist(),
            },
            "result": {
                "TI": np.zeros((6, 8, len(ns.frequencies_ghz)), dtype=np.float32),
                "TV": np.zeros((6, 8, len(ns.frequencies_ghz)), dtype=np.float32),
            },
            "outputs": {
                "output_dir": str(tmp_path),
                "save_outputs": False,
                "write_preview": False,
                "h5_path": None,
                "preview_png": None,
                "fits_paths": [],
            },
        }

    monkeypatch.setattr(render_mw_workflow, "run", _fake_run)

    result = render_mw_maps(
        MWRenderOptions(
            model_path=tmp_path / "model.chr.h5",
            model_format="h5",
            freqlist_ghz=[5.7, 8.1, 11.2],
            plasma=CoronalPlasmaParameters(
                tbase=3.0e6,
                nbase=2.5e8,
                q0=0.025,
                a=0.6,
                b=1.9,
                mode=5,
                shtable=custom_shtable,
            ),
            save_outputs=False,
            write_preview=False,
            verbose=True,
        )
    )

    assert captured["verbose"] is True
    assert captured["ns"].frequencies_ghz == [5.7, 8.1, 11.2]
    assert np.isclose(captured["ns"].tbase, 3.0e6)
    assert np.isclose(captured["ns"].nbase, 2.5e8)
    assert np.isclose(captured["ns"].q0, 0.025)
    assert np.isclose(captured["ns"].a, 0.6)
    assert np.isclose(captured["ns"].b, 1.9)
    assert captured["ns"].corona_mode == 5
    assert np.array_equal(np.asarray(captured["ns"].shtable), custom_shtable)
    assert result.freqlist_ghz == [5.7, 8.1, 11.2]
    assert np.isclose(result.plasma.tbase, 3.0e6)
    assert np.isclose(result.plasma.nbase, 2.5e8)
    assert np.isclose(result.plasma.q0, 0.025)
    assert np.isclose(result.plasma.a, 0.6)
    assert np.isclose(result.plasma.b, 1.9)
    assert result.plasma.mode == 5
    assert np.array_equal(np.asarray(result.plasma.shtable), custom_shtable)
