"""Sphinx configuration for gximagecomputing."""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

project = "pyGXrender / gximagecomputing"
author = "suncast-org"
copyright = "2026, suncast-org"
release = "0.0.2.1"

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent
SRC_DIR = REPO_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

# Avoid SunPy/Matplotlib config/cache permission issues during autodoc imports.
_tmp_root = Path(tempfile.gettempdir()) / "gximagecomputing_sphinx"
for _name in ("sunpy_cfg", "mpl_cfg", "xdg_cache"):
    (_tmp_root / _name).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("SUNPY_CONFIGDIR", str(_tmp_root / "sunpy_cfg"))
os.environ.setdefault("MPLCONFIGDIR", str(_tmp_root / "mpl_cfg"))
os.environ.setdefault("XDG_CACHE_HOME", str(_tmp_root / "xdg_cache"))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

autosummary_generate = True
autoclass_content = "both"
autodoc_member_order = "bysource"
autodoc_typehints = "description"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Keep docs builds focused on actionable warnings from our sources.
suppress_warnings = [
    "toc.excluded",
]
