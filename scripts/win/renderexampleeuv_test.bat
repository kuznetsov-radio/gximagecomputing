@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..") do set "REPO_ROOT=%%~fI"

set "OUTDIR=%TEMP%\gximagecomputing_validation_groundtruth"
set "OUTNAME=test.chr.h5_py_euv_maps.h5"
set "MODEL_NAME=test.chr.h5"
set "EBTEL_NAME=ebtel.sav"
set "INSTRUMENT=aia"

set "PYTHONPATH=src"

for /f "usebackq delims=" %%A in (`python -m gxrender.utils.test_data model "%MODEL_NAME%"`) do set "MODEL=%%A"
if not defined MODEL (
  echo Failed to resolve model test data.
  exit /b 1
)

for /f "usebackq delims=" %%A in (`python -m gxrender.utils.test_data response "%INSTRUMENT%"`) do set "RESPONSE=%%A"
if not defined RESPONSE (
  echo Failed to resolve response test data.
  exit /b 1
)

for /f "usebackq delims=" %%A in (`python -m gxrender.utils.test_data ebtel "%EBTEL_NAME%"`) do set "EBTEL=%%A"
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

python src\gxrender\workflows\render_euv.py ^
  --model-path "%MODEL%" ^
  --instrument "%INSTRUMENT%" ^
  --response-sav "%RESPONSE%" ^
  --ebtel-path "%EBTEL%" ^
  --output-dir "%OUTDIR%" ^
  --output-name "%OUTNAME%"
if errorlevel 1 exit /b %errorlevel%

echo You may use gxrender-map-view "%OUTDIR%\%OUTNAME%" to visualize the results
