@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..") do set "REPO_ROOT=%%~fI"

set "OUTDIR=%TEMP%\gximagecomputing_validation_groundtruth"
set "OUTNAME=test.chr.h5_py_mw_maps.h5"

set "PYTHONPATH=src"

for /f "usebackq delims=" %%A in (`python -m gxrender.utils.test_data model test.chr.h5`) do set "MODEL=%%A"
if not defined MODEL (
  echo Failed to resolve model test data.
  exit /b 1
)

for /f "usebackq delims=" %%A in (`python -m gxrender.utils.test_data ebtel ebtel.sav`) do set "EBTEL=%%A"
if not defined EBTEL (
  echo Failed to resolve EBTEL test data.
  exit /b 1
)

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
if not exist "%TEMP%\gximagecomputing_sunpy" mkdir "%TEMP%\gximagecomputing_sunpy"
if not exist "%TEMP%\gximagecomputing_mpl" mkdir "%TEMP%\gximagecomputing_mpl"

cd /d "%REPO_ROOT%"
set "SUNPY_CONFIGDIR=%TEMP%\gximagecomputing_sunpy"
set "MPLCONFIGDIR=%TEMP%\gximagecomputing_mpl"

python src\gxrender\workflows\render_mw.py ^
  --model-path "%MODEL%" ^
  --ebtel-path "%EBTEL%" ^
  --output-dir "%OUTDIR%" ^
  --output-name "%OUTNAME%"
if errorlevel 1 exit /b %errorlevel%

echo You may use gxrender-map-view "%OUTDIR%\%OUTNAME%" to visualize the results
