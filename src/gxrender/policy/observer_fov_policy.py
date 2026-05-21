from __future__ import annotations

from gxrender.policy.contracts import ObserverFovRequest, ObserverFovResolution
from gxrender.workflows import _render_common


def resolve_observer_fov_policy(request: ObserverFovRequest) -> ObserverFovResolution:
    """Phase-1 wrapper around existing observer/FOV logic (no behavior changes)."""
    (
        model,
        model_dt,
        model_metadata,
        observer_geometry,
        center_source,
        xc_auto,
        yc_auto,
        model_w_arcsec,
        model_h_arcsec,
        applied_overrides,
    ) = _render_common._load_model_and_fov_impl(
        model_path=request.model_path,
        loader=request.loader,
        cli_args=request.cli_args,
        prefer_execute_center=request.prefer_execute_center,
    )

    return ObserverFovResolution(
        model=model,
        model_dt=model_dt,
        model_metadata=dict(model_metadata),
        observer_geometry=observer_geometry,
        center_source=str(center_source),
        xc_auto=float(xc_auto),
        yc_auto=float(yc_auto),
        model_w_arcsec=float(model_w_arcsec),
        model_h_arcsec=float(model_h_arcsec),
        applied_overrides=dict(applied_overrides),
        decisions=[],
    )
