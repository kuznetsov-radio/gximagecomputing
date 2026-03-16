# Utilities for gxrender Python workflows.

__all__ = [
    "find_ebtel_file",
    "find_model_file",
    "find_response_file",
    "test_data_setup_hint",
    "try_find_ebtel_file",
    "try_find_model_file",
    "try_find_response_file",
]


def __getattr__(name):
    if name in __all__:
        from . import test_data as _test_data

        return getattr(_test_data, name)
    raise AttributeError(f"module 'gxrender.utils' has no attribute {name!r}")
