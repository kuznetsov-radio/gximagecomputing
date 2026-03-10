"""Compatibility wrapper for rendered-map viewer.

Canonical module path is pyGXrender.viz.render_map_view.
Implementation currently lives in pyGXrender.utils.render_map_view.
"""

from pyGXrender.utils.render_map_view import main, parse_args, run_viewer

__all__ = ["main", "parse_args", "run_viewer"]
