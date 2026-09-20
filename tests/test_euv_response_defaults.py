from __future__ import annotations

from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from gxrender import euv as gx_euv
from gxrender import sdk as gx_sdk
from gxrender.euv import EUVResponseMeta
from gxrender.geometry.observer_geometry import ResolvedObserverGeometry
from gxrender.workflows import render_euv


def _base_args(**overrides):
    values = {
        "instrument": "AIA",
        "channels": ["171"],
        "response": None,
        "response_dt": None,
        "response_meta": None,
        "response_sav": None,
    }
    values.update(overrides)
    return Namespace(**values)


def _observer_geometry(*, observer_name: str, observer_source: str) -> ResolvedObserverGeometry:
    return ResolvedObserverGeometry(
        observer_name=observer_name,
        l0_deg=0.0,
        b0_deg=0.0,
        dsun_cm=1.0,
        render_lonc_deg=0.0,
        render_b0_deg=0.0,
        render_dsun_cm=1.0,
        observer_source=observer_source,
        warnings=(),
        rsun_cm=None,
        rsun_arcsec=None,
    )


def test_default_instrument_infers_stereo_a_from_observer_metadata() -> None:
    args = _base_args(instrument=None)

    render_euv._apply_default_response_selection(
        args,
        observer_geometry=_observer_geometry(observer_name="stereo-a", observer_source="saved_observer_metadata"),
    )

    assert args.instrument == "STEREO-A"


def test_default_instrument_rejects_ambiguous_solar_orbiter_observer() -> None:
    args = _base_args(instrument=None)

    render_euv._apply_default_response_selection(
        args,
        observer_geometry=_observer_geometry(observer_name="solar orbiter", observer_source="saved_observer_metadata"),
    )

    assert args.instrument == "SOLO-FSI"


def test_default_instrument_falls_back_to_aia_only_without_observer_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    args = _base_args(instrument=None)
    messages: list[str] = []
    monkeypatch.setattr(render_euv, "_warn_example_default", lambda message: messages.append(message))

    render_euv._apply_default_response_selection(
        args,
        observer_geometry=_observer_geometry(observer_name="earth", observer_source="default_earth"),
    )

    assert args.instrument == "AIA"
    assert any("AIA was assumed" in message for message in messages)


def test_explicit_instrument_override_is_preserved() -> None:
    args = _base_args(instrument="SOLO-HRI")

    render_euv._apply_default_response_selection(
        args,
        observer_geometry=_observer_geometry(observer_name="solar orbiter", observer_source="saved_observer_metadata"),
    )

    assert args.instrument == "SOLO-HRI"


def test_resolve_response_inputs_prefers_python_native_aia_provider(monkeypatch: pytest.MonkeyPatch) -> None:
    payload = np.zeros(1, dtype=[("ds", np.float64), ("NT", np.int32), ("Nchannels", np.int32)])
    response_dt = payload.dtype
    response_meta = EUVResponseMeta(
        instrument="AIA",
        channels=["A171"],
        source="pyeuvtools-export",
        mode="pyeuvtools:evenorm_chiantifix",
    )

    received = {}

    def fake_builder(**kwargs):
        received.update(kwargs)
        return payload, response_dt, response_meta

    monkeypatch.setattr(
        render_euv,
        "build_default_aia_euv_response",
        fake_builder,
    )
    monkeypatch.setattr(
        render_euv,
        "_resolve_default_response_sav",
        lambda instrument: pytest.fail(f"unexpected SAV lookup for {instrument}"),
    )

    response, resolved_dt, resolved_meta = render_euv._resolve_response_inputs(
        _base_args(),
        obs_time_iso="2025-11-26T15:34:31",
    )

    assert response is payload
    assert resolved_dt is response_dt
    assert resolved_meta is response_meta
    assert received["correction_state"] == "evenorm_chiantifix"


def test_default_aia_response_requests_idl_default_corrections(monkeypatch: pytest.MonkeyPatch) -> None:
    payload = np.zeros(1, dtype=[("ds", np.float64), ("NT", np.int32), ("Nchannels", np.int32)])
    captured = {}

    class FakeExport:
        channels = ("94", "131", "171", "193", "211", "304", "335")
        emissivity_wavelength = np.asarray([10.0])
        emissivity_logte = np.asarray([6.0])
        emissivity = np.asarray([[1.0]])

    def fake_payload_builder(**kwargs):
        captured.update(kwargs)
        return payload, payload.dtype, {
            "instrument": "AIA",
            "channels": tuple(f"A{channel}" for channel in FakeExport.channels),
        }

    monkeypatch.setattr(
        gx_euv,
        "_require_pyeuvtools_aia_bridge",
        lambda: (fake_payload_builder, lambda path: FakeExport(), lambda: Path("hybrid.sav")),
    )

    _response, _response_dt, metadata = gx_euv.build_default_aia_euv_response(
        obstime="2012-07-12T04:46:25.800",
    )

    assert captured["channels"] == list(FakeExport.channels)
    assert captured["include_eve_correction"] is True
    assert captured["include_chiantifix"] is True
    assert metadata.mode == "pyeuvtools:evenorm_chiantifix"


def test_resolve_response_inputs_honors_explicit_response_sav(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    sav_path = tmp_path / "resp_aia.sav"
    loaded = object()
    loaded_dt = object()
    loaded_meta = EUVResponseMeta(instrument="AIA", channels=["A171"], source=str(sav_path), mode="sav")

    monkeypatch.setattr(
        render_euv,
        "build_default_aia_euv_response",
        lambda **kwargs: pytest.fail(f"unexpected Python-native provider call: {kwargs}"),
    )
    monkeypatch.setattr(
        render_euv,
        "load_euv_response_sav",
        lambda path: (loaded, loaded_dt, loaded_meta) if path == str(sav_path) else pytest.fail(path),
    )

    response, response_dt, response_meta = render_euv._resolve_response_inputs(
        _base_args(response_sav=sav_path),
        obs_time_iso="2025-11-26T15:34:31",
    )

    assert response is loaded
    assert response_dt is loaded_dt
    assert response_meta is loaded_meta


def test_sdk_result_exposes_response_source_metadata(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    workflow_result = {
        "library_path": "/tmp/libgxeuv.dylib",
        "model_path": str(tmp_path / "model.h5"),
        "model_format": "h5",
        "ebtel_path": "",
        "observer_overrides_applied": {},
        "center_source": "default",
        "geometry": {
            "xc_arcsec": 0.0,
            "yc_arcsec": 0.0,
            "dx_arcsec": 2.0,
            "dy_arcsec": 2.0,
            "nx": 4,
            "ny": 3,
            "fov_x_arcsec": 8.0,
            "fov_y_arcsec": 6.0,
        },
        "obs_time_iso": "2025-11-26T15:34:31",
        "plasma": {
            "tbase_k": 1.0e6,
            "nbase_cm3": 1.0e8,
            "q0": 0.0217,
            "a": 0.3,
            "b": 2.7,
            "mode": 0,
            "selective_heating": False,
            "shtable": None,
        },
        "response": {
            "instrument": "AIA",
            "channels": ["A171"],
            "source": "pyeuvtools-export",
            "mode": "pyeuvtools:evenorm_chiantifix",
        },
        "result": {
            "flux_corona": np.zeros((1, 3, 4), dtype=np.float64),
            "flux_tr": np.zeros((1, 3, 4), dtype=np.float64),
        },
        "outputs": {
            "output_dir": str(tmp_path),
            "save_outputs": False,
            "write_preview": False,
            "h5_path": None,
            "preview_png": None,
        },
    }

    monkeypatch.setattr(gx_sdk._render_euv_workflow, "run", lambda ns, verbose=False: workflow_result)

    result = gx_sdk.render_euv_maps(
        gx_sdk.EUVRenderOptions(
            model_path=tmp_path / "model.h5",
            model_format="h5",
            instrument="AIA",
            channels=["171"],
            save_outputs=False,
            write_preview=False,
        )
    )

    assert result.response.instrument == "AIA"
    assert result.response.channels == ["A171"]
    assert result.response.source == "pyeuvtools-export"
    assert result.response.mode == "pyeuvtools:evenorm_chiantifix"


def test_sdk_forwards_named_observer_to_euv_workflow(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    observed = {}

    def fake_run(ns, verbose=False):
        observed["observer"] = getattr(ns, "observer", None)
        observed["parallel"] = getattr(ns, "parallel", None)
        observed["exact"] = getattr(ns, "exact", None)
        observed["projection_threads"] = getattr(ns, "projection_threads", None)
        return {
            "library_path": "/tmp/libgxeuv.dylib",
            "model_path": str(tmp_path / "model.h5"),
            "model_format": "h5",
            "ebtel_path": "",
            "observer_overrides_applied": {},
            "center_source": "default",
            "geometry": {
                "xc_arcsec": 0.0,
                "yc_arcsec": 0.0,
                "dx_arcsec": 2.0,
                "dy_arcsec": 2.0,
                "nx": 4,
                "ny": 3,
                "fov_x_arcsec": 8.0,
                "fov_y_arcsec": 6.0,
            },
            "obs_time_iso": "2025-11-26T15:34:31",
            "plasma": {
                "tbase_k": 1.0e6,
                "nbase_cm3": 1.0e8,
                "q0": 0.0217,
                "a": 0.3,
                "b": 2.7,
                "mode": 0,
                "selective_heating": False,
                "shtable": None,
            },
            "response": {
                "instrument": "SOLO-FSI",
                "channels": ["174"],
                "source": "resp_solo-fsi.sav",
                "mode": "sav",
            },
            "result": {
                "flux_corona": np.zeros((1, 3, 4), dtype=np.float64),
                "flux_tr": np.zeros((1, 3, 4), dtype=np.float64),
            },
            "outputs": {
                "output_dir": str(tmp_path),
                "save_outputs": False,
                "write_preview": False,
                "h5_path": None,
                "preview_png": None,
            },
        }

    monkeypatch.setattr(gx_sdk._render_euv_workflow, "run", fake_run)

    gx_sdk.render_euv_maps(
        gx_sdk.EUVRenderOptions(
            model_path=tmp_path / "model.h5",
            model_format="h5",
            observer_name="solo",
            parallel=True,
            exact=False,
            projection_threads=8,
            save_outputs=False,
            write_preview=False,
        )
    )

    assert observed["observer"] == "solo"
    assert observed["parallel"] is True
    assert observed["exact"] is False
    assert observed["projection_threads"] == 8
