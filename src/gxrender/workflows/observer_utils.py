"""Compatibility wrappers for legacy observer workflow imports.

This module historically carried a second copy of the observer-geometry logic.
Keep the import path stable, but delegate to the canonical implementation in
``gxrender.geometry.observer_geometry`` so the workflow layer no longer owns a
separate geometry stack.
"""

from gxrender.geometry.observer_geometry import (
    MetadataObserverState,
    ResolvedObserverGeometry,
    build_ephemeris_from_pb0r,
    build_observer_coordinate,
    build_pb0r_from_ephemeris,
    compute_inscribing_fov,
    compute_inscribing_fov_box,
    compute_projected_fov_for_observer,
    compute_sunpy_wcs_header,
    model_time_from_model,
    normalize_observer_name,
    observer_summary,
    resolve_observer_geometry,
    resolve_simbox_from_observer_and_model,
    should_use_saved_observer_fov,
)

__all__ = [
    "MetadataObserverState",
    "ResolvedObserverGeometry",
    "build_ephemeris_from_pb0r",
    "build_observer_coordinate",
    "build_pb0r_from_ephemeris",
    "compute_inscribing_fov",
    "compute_inscribing_fov_box",
    "compute_projected_fov_for_observer",
    "compute_sunpy_wcs_header",
    "model_time_from_model",
    "normalize_observer_name",
    "observer_summary",
    "resolve_observer_geometry",
    "resolve_simbox_from_observer_and_model",
    "should_use_saved_observer_fov",
]
