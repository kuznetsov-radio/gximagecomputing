from __future__ import annotations

from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from gxrender.utils.test_data import (
    try_find_model_loader_parity_files,
    try_find_response_file,
)
from gxrender.workflows import _render_common, render_euv, render_mw


pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:Current Python MW workflow uses projection=0.*:UserWarning"
    ),
    pytest.mark.filterwarnings(
        "ignore:Current Python EUV workflow uses projection flags off.*:UserWarning"
    ),
]


@pytest.fixture(scope="module")
def h5_model_path() -> Path:
    pair = try_find_model_loader_parity_files()
    if pair is None:
        pytest.skip("loader parity fixtures are unavailable")
    return pair[1]


def _mw_args(model_path: Path) -> Namespace:
    return Namespace(
        model_path=model_path,
        model_format="auto",
        ebtel_path="",
        output_dir=Path("/tmp"),
        output_name="mw_policy_parity.h5",
        output_format="h5",
        omp_threads=1,
        xc=None,
        yc=None,
        dsun_cm=None,
        lonc_deg=None,
        b0sun_deg=None,
        observer=None,
        dx=6.0,
        dy=6.0,
        pixel_scale_arcsec=None,
        nx=32,
        ny=32,
        xrange=None,
        yrange=None,
        auto_fov=False,
        use_saved_fov=False,
        frequencies_ghz=[5.8],
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


def _euv_args(model_path: Path, response_sav: Path | None) -> Namespace:
    if response_sav is None:
        pytest.skip("AIA response fixture is unavailable for EUV parity integration test")
    return Namespace(
        model_path=model_path,
        model_format="auto",
        ebtel_path="",
        output_dir=Path("/tmp"),
        output_name="euv_policy_parity.h5",
        channels=["171"],
        instrument=None,
        response_sav=response_sav,
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
        dx=6.0,
        dy=6.0,
        pixel_scale_arcsec=None,
        nx=32,
        ny=32,
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


def test_mw_full_output_parity_wrapper_vs_legacy(h5_model_path: Path, monkeypatch):
    args = _mw_args(h5_model_path)

    routed = render_mw.run(args, verbose=False)

    monkeypatch.setattr(
        _render_common,
        "load_model_and_fov",
        lambda model_path, loader, cli_args, prefer_execute_center=True: _render_common._load_model_and_fov_impl(
            model_path,
            loader,
            cli_args,
            prefer_execute_center=prefer_execute_center,
        ),
    )
    legacy = render_mw.run(args, verbose=False)

    assert routed["geometry"] == legacy["geometry"]
    assert routed["center_source"] == legacy["center_source"]
    assert routed["freqlist_ghz"] == legacy["freqlist_ghz"]
    assert np.allclose(routed["result"]["TI"], legacy["result"]["TI"], rtol=0.0, atol=0.0)
    assert np.allclose(routed["result"]["TV"], legacy["result"]["TV"], rtol=0.0, atol=0.0)


def test_euv_full_output_parity_wrapper_vs_legacy(h5_model_path: Path, monkeypatch):
    args = _euv_args(h5_model_path, try_find_response_file("aia"))

    routed = render_euv.run(args, verbose=False)

    def _legacy_response_wrapper(request):
        response, response_dt, response_meta = render_euv._resolve_response_inputs(
            request.args,
            obs_time_iso=request.obs_time_iso,
        )
        return SimpleNamespace(
            response=response,
            response_dt=response_dt,
            response_meta=response_meta,
        )

    monkeypatch.setattr(
        render_euv,
        "resolve_euv_response",
        _legacy_response_wrapper,
    )
    legacy = render_euv.run(args, verbose=False)

    assert routed["geometry"] == legacy["geometry"]
    assert routed["center_source"] == legacy["center_source"]
    assert routed["response"]["instrument"] == legacy["response"]["instrument"]
    assert routed["response"]["channels"] == legacy["response"]["channels"]
    assert np.allclose(routed["result"]["flux_corona"], legacy["result"]["flux_corona"], rtol=0.0, atol=0.0)
    assert np.allclose(routed["result"]["flux_tr"], legacy["result"]["flux_tr"], rtol=0.0, atol=0.0)
