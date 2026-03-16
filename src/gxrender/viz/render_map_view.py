"""Compatibility wrapper for rendered-map viewer.

Canonical module path is gxrender.viz.render_map_view.
Implementation currently lives in gxrender.utils.render_map_view.
"""

from gxrender.utils.render_map_view import main, parse_args, run_viewer

__all__ = ["main", "parse_args", "run_viewer"]
