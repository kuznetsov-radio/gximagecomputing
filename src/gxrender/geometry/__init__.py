"""Shared observer and geometry helpers for gxrender."""

from .observer_geometry import (
    ResolvedObserverGeometry,
    build_ephemeris_from_pb0r,
    build_observer_coordinate,
    build_pb0r_from_ephemeris,
    compute_inscribing_fov,
    compute_inscribing_fov_box,
    compute_projected_fov_for_observer,
    compute_sunpy_wcs_header,
    normalize_observer_name,
    observer_summary,
    resolve_observer_geometry,
    resolve_simbox_from_observer_and_model,
)

__all__ = [
    "ResolvedObserverGeometry",
    "build_ephemeris_from_pb0r",
    "build_observer_coordinate",
    "build_pb0r_from_ephemeris",
    "compute_inscribing_fov",
    "compute_inscribing_fov_box",
    "compute_projected_fov_for_observer",
    "compute_sunpy_wcs_header",
    "normalize_observer_name",
    "observer_summary",
    "resolve_observer_geometry",
    "resolve_simbox_from_observer_and_model",
]
