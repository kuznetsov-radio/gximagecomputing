from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import pytest

from gxrender.policy.contracts import EUVResponseRequest, ObserverFovRequest
from gxrender.policy.euv_response_policy import apply_default_response_selection, resolve_euv_response
from gxrender.policy.observer_fov_policy import resolve_observer_fov_policy
from gxrender.policy.response_providers import (
    clear_registered_euv_response_providers,
    register_euv_response_provider,
)
from gxrender.policy.script_policy import evaluate_script_policy
from gxrender.workflows import _render_common, render_euv


pytestmark = [
    pytest.mark.filterwarnings("ignore:Example EUV default applied.*:UserWarning"),
    pytest.mark.filterwarnings(
        "ignore:Current Python EUV workflow uses projection flags off.*:UserWarning"
    ),
]


def test_observer_fov_wrapper_parity_with_legacy(monkeypatch):
    observer_geometry = SimpleNamespace(observer_source="default_earth")
    expected_tuple = (
        {"model": "payload"},
        "model_dt",
        {"obs_time": "2020-01-01T00:00:00"},
        observer_geometry,
        "saved_fov",
        -10.0,
        20.0,
        300.0,
        400.0,
        {"DSun_cm": 1.5e13},
    )

    def _fake_load_model_and_fov_impl(*args, **kwargs):
        return expected_tuple

    monkeypatch.setattr(_render_common, "_load_model_and_fov_impl", _fake_load_model_and_fov_impl)

    request = ObserverFovRequest(
        model_path=Path("dummy.h5"),
        loader="h5",
        cli_args=Namespace(),
        prefer_execute_center=True,
    )

    legacy = _render_common._load_model_and_fov_impl(
        model_path=request.model_path,
        loader=request.loader,
        cli_args=request.cli_args,
        prefer_execute_center=request.prefer_execute_center,
    )
    wrapped = resolve_observer_fov_policy(request)

    assert wrapped.model == legacy[0]
    assert wrapped.model_dt == legacy[1]
    assert wrapped.model_metadata == legacy[2]
    assert wrapped.observer_geometry == legacy[3]
    assert wrapped.center_source == legacy[4]
    assert wrapped.xc_auto == legacy[5]
    assert wrapped.yc_auto == legacy[6]
    assert wrapped.model_w_arcsec == legacy[7]
    assert wrapped.model_h_arcsec == legacy[8]
    assert wrapped.applied_overrides == legacy[9]


def test_load_model_and_fov_routes_through_policy_wrapper(monkeypatch):
    observer_geometry = SimpleNamespace(observer_source="default_earth")

    def _fake_resolve_observer_fov_policy(request):
        return SimpleNamespace(
            model={"model": "payload"},
            model_dt="model_dt",
            model_metadata={"obs_time": "2020-01-01T00:00:00"},
            observer_geometry=observer_geometry,
            center_source="saved_fov",
            xc_auto=-10.0,
            yc_auto=20.0,
            model_w_arcsec=300.0,
            model_h_arcsec=400.0,
            applied_overrides={"DSun_cm": 1.5e13},
        )

    monkeypatch.setattr(
        "gxrender.policy.observer_fov_policy.resolve_observer_fov_policy",
        _fake_resolve_observer_fov_policy,
    )

    out = _render_common.load_model_and_fov(
        model_path=Path("dummy.h5"),
        loader="h5",
        cli_args=Namespace(),
        prefer_execute_center=True,
    )

    assert out[0] == {"model": "payload"}
    assert out[4] == "saved_fov"
    assert out[5] == -10.0
    assert out[6] == 20.0


def test_apply_default_response_selection_wrapper_parity():
    observer_geometry = SimpleNamespace(observer_name="earth", observer_source="default_earth")

    args_legacy = Namespace(response=None, response_sav=None, instrument=None)
    args_wrapped = Namespace(response=None, response_sav=None, instrument=None)

    render_euv._apply_default_response_selection(args_legacy, observer_geometry=observer_geometry)
    apply_default_response_selection(args_wrapped, observer_geometry=observer_geometry)

    assert args_wrapped.instrument == args_legacy.instrument


def test_resolve_euv_response_wrapper_parity_for_prebuilt_response():
    response = object()
    response_dt = object()
    response_meta = SimpleNamespace(source="prebuilt", mode="object", instrument="AIA", channels=["171"])

    args = Namespace(
        response=response,
        response_dt=response_dt,
        response_meta=response_meta,
        response_sav=None,
        instrument=None,
        channels=None,
    )

    legacy = render_euv._resolve_response_inputs(args, obs_time_iso="2020-01-01T00:00:00")
    wrapped = resolve_euv_response(EUVResponseRequest(args=args, obs_time_iso="2020-01-01T00:00:00"))

    assert wrapped.response is legacy[0]
    assert wrapped.response_dt is legacy[1]
    assert wrapped.response_meta is legacy[2]
    assert wrapped.source == "prebuilt"
    assert wrapped.mode == "object"


def test_render_euv_run_uses_policy_wrappers(monkeypatch):
    args = Namespace(
        model_path=Path("dummy.h5"),
        model_format="h5",
        ebtel_path="",
        output_dir=Path("/tmp"),
        output_name="dummy.h5",
        channels=None,
        instrument=None,
        response_sav=None,
        response=None,
        response_dt=None,
        response_meta=None,
        omp_threads=1,
        xc=None,
        yc=None,
        dsun_cm=None,
        lonc_deg=None,
        b0sun_deg=None,
        observer=None,
        dx=2.0,
        dy=2.0,
        pixel_scale_arcsec=None,
        nx=16,
        ny=16,
        xrange=None,
        yrange=None,
        auto_fov=False,
        use_saved_fov=False,
        tbase=1e6,
        nbase=1e8,
        q0=0.0217,
        a=0.3,
        b=2.7,
        shtable_path=None,
        selective_heating=False,
        corona_mode=0,
        force_isothermal=False,
        interpol_b=False,
        analytical_nt=False,
        save_outputs=False,
        write_preview=False,
    )

    observer_geometry = SimpleNamespace(
        observer_name="earth",
        observer_source="default_earth",
        warnings=[],
        l0_deg=0.0,
        b0_deg=0.0,
        dsun_cm=1.5e13,
        rsun_cm=6.96e10,
        rsun_arcsec=960.0,
    )
    common = SimpleNamespace(
        model_path=Path("dummy.h5"),
        loader="h5",
        model={"obstime": [0.0]},
        model_dt="model_dt",
        ebtel_path="",
        ebtel_c="ebtel_c",
        ebtel_dt="ebtel_dt",
        center_source="saved_fov",
        xc=0.0,
        yc=0.0,
        dx=2.0,
        dy=2.0,
        nx=16,
        ny=16,
        fov_x=32.0,
        fov_y=32.0,
        observer_geometry=observer_geometry,
        observer_overrides_applied={},
    )

    calls = {"default": 0, "resolve": 0}

    monkeypatch.setattr(render_euv, "prepare_common_inputs", lambda *_args, **_kwargs: common)
    monkeypatch.setattr(
        render_euv,
        "resolve_plasma_parameters",
        lambda *_args, **_kwargs: SimpleNamespace(
            tbase=1e6,
            nbase=1e8,
            q0=0.0217,
            a=0.3,
            b=2.7,
            mode=0,
            shtable=None,
            selective_heating=False,
        ),
    )
    monkeypatch.setattr(render_euv, "observer_summary", lambda *_args, **_kwargs: "earth")
    monkeypatch.setattr(render_euv, "model_obstime_iso", lambda *_args, **_kwargs: "2020-01-01T00:00:00")
    monkeypatch.setattr(render_euv, "compute_sunpy_wcs_header", lambda **_kwargs: {"DATE-OBS": "2020-01-01T00:00:00"})

    class _FakeGX:
        libname = "fake_lib"

        def synth_euv(self, **_kwargs):
            import numpy as np

            return {
                "flux_corona": np.zeros((16, 16, 1), dtype=float),
                "flux_tr": np.zeros((16, 16, 1), dtype=float),
            }

    monkeypatch.setattr(render_euv, "GXEUVImageComputing", _FakeGX)

    def _fake_apply_default_response_selection(args_value, *, observer_geometry):
        calls["default"] += 1
        args_value.instrument = "AIA"
        return args_value

    def _fake_resolve_euv_response(request):
        calls["resolve"] += 1
        meta = SimpleNamespace(instrument="AIA", channels=["171"], source="prebuilt", mode="object")
        return SimpleNamespace(response="response", response_dt="response_dt", response_meta=meta)

    monkeypatch.setattr(render_euv, "apply_default_response_selection", _fake_apply_default_response_selection)
    monkeypatch.setattr(render_euv, "resolve_euv_response", _fake_resolve_euv_response)

    out = render_euv.run(args, verbose=False)

    assert calls["default"] == 1
    assert calls["resolve"] == 1
    assert out["response"]["instrument"] == "AIA"


def test_provider_registry_uses_registered_provider_without_flag(monkeypatch):
    class _Provider:
        name = "mock-provider"
        priority = 1

        def supports(self, instrument, obs_time_iso):
            return instrument == "stereo-a"

        def build(self, instrument, obs_time_iso, channels):
            meta = SimpleNamespace(source="mock-provider", mode="python_provider_time_dependent")
            return "response", "response_dt", meta

    args = Namespace(
        response=None,
        response_dt=None,
        response_meta=None,
        response_sav=None,
        instrument="STEREO-A",
        channels=["171"],
    )

    clear_registered_euv_response_providers()
    register_euv_response_provider(_Provider())

    wrapped = resolve_euv_response(EUVResponseRequest(args=args, obs_time_iso="2020-01-01T00:00:00"))

    assert wrapped.response == "response"
    assert wrapped.response_dt == "response_dt"
    assert wrapped.source == "mock-provider"
    assert wrapped.mode == "python_provider_time_dependent"
    assert wrapped.decisions == ["provider:mock-provider"]

    clear_registered_euv_response_providers()


def test_provider_registry_falls_back_to_legacy_for_explicit_prebuilt_response(monkeypatch):
    args = Namespace(
        response=object(),
        response_dt=object(),
        response_meta=SimpleNamespace(source="prebuilt", mode="object", instrument="AIA", channels=["171"]),
        response_sav=None,
        instrument="STEREO-A",
        channels=["171"],
    )

    wrapped = resolve_euv_response(EUVResponseRequest(args=args, obs_time_iso="2020-01-01T00:00:00"))

    assert wrapped.source == "prebuilt"
    assert wrapped.mode == "object"


def test_script_policy_allows_implicit_non_earth_when_provider_exists(monkeypatch):
    class _Provider:
        name = "mock-provider"
        priority = 1

        def supports(self, instrument, obs_time_iso):
            return instrument == "stereo-a"

        def build(self, instrument, obs_time_iso, channels):
            raise AssertionError("build should not be called during policy evaluation")

    clear_registered_euv_response_providers()
    register_euv_response_provider(_Provider())

    monkeypatch.setattr(
        "gxrender.policy.script_policy._load_model_metadata",
        lambda _path: ({"observer_name": "earth", "obs_time": "2020-01-01T00:00:00"}, None),
    )

    result = evaluate_script_policy(
        model_path=Path("dummy.h5"),
        instrument="",
        observer="stereo-a",
        response_sav="",
    )

    assert result.normalized_observer == "stereo-a"
    assert result.implicit_response_mode == "provider_time_dependent"

    clear_registered_euv_response_providers()
