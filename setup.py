#!/usr/bin/env python3

import platform
from setuptools import setup, Extension
from pathlib import Path
import os
import subprocess

common_compile_flags = ["-std=c++11", "-O3", "-fPIC"]

compile_flags = {
    "Linux": ["-DLINUX", *common_compile_flags],
    "Windows": ["-DWINDOWS"],
    "Darwin": ["-DMACOS", *common_compile_flags, "-Xpreprocessor", "-fopenmp"],
}

link_flags = {
    "Linux": ["-fopenmp"],
    "Windows": ["-fopenmp"],
    "Darwin": ["-lomp"],
}

current_os = platform.system()
print("building for", current_os)


def _env_flag(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


# Default policy:
# - Windows: skip local native compilation unless explicitly requested.
# - Linux/macOS: keep compiling by default.
build_native = _env_flag("PYGXRENDER_BUILD_NATIVE", default=(current_os != "Windows"))

if current_os in link_flags:
    extra_link = link_flags[current_os]
    extra_compile = compile_flags[current_os]
else:
    extra_link = []
    extra_compile = []

if current_os == "Darwin":
    brew_prefix = os.environ.get("HOMEBREW_PREFIX")
    if not brew_prefix:
        try:
            brew_prefix = subprocess.check_output(
                ["brew", "--prefix"], text=True
            ).strip()
        except Exception:
            brew_prefix = None
    if brew_prefix:
        extra_compile += [
            f"-I{brew_prefix}/include",
            f"-I{brew_prefix}/opt/libomp/include",
        ]
        extra_link += [f"-L{brew_prefix}/lib", f"-L{brew_prefix}/opt/libomp/lib"]

source_dir = Path("./source")
source_files = [str(source_dir / x.name) for x in sorted(source_dir.glob("*.cpp"))]

if current_os != "Windows":
    source_files.remove(str(source_dir / "dllmain.cpp"))

ext_modules = []
if build_native:
    source_files.append("pyinit.cpp")
    # Define the extension module
    render_grff_module = Extension(
        'gxrender.RenderGRFF',
        sources=source_files,
        extra_compile_args=extra_compile,
        extra_link_args=extra_link,
        include_dirs=[source_dir],
        export_symbols=[],
        language="c++",
    )
    ext_modules = [render_grff_module]
else:
    print(
        "Skipping native extension build "
        "(set PYGXRENDER_BUILD_NATIVE=1 to force local compilation)."
    )

setup(
    name="pyGXrender",
    ext_modules=ext_modules,
    include_package_data=True
)
