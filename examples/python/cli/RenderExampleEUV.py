#!/usr/bin/env python3

"""Compatibility CLI wrapper for the EUV rendering workflow.

This script remains in examples/ for user familiarity, but the implementation
now lives in the installable package module:
    gximagecomputing.workflows.render_euv
"""

from gximagecomputing.workflows.render_euv import main


if __name__ == "__main__":
    main()
