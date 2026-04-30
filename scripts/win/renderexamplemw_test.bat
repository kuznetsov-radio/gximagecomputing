@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..") do set "REPO_ROOT=%%~fI"

set "OUTDIR=%TEMP%\gximagecomputing_validation_groundtruth"
if not defined OUTNAME set "OUTNAME="
if not defined MODEL_NAME set "MODEL_NAME="
if not defined EBTEL_NAME set "EBTEL_NAME=ebtel.sav"
if defined PYTHON_BIN (
  set "PYTHON_CMD=%PYTHON_BIN%"
) else (
  set "PYTHON_CMD=python"
)

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
if not exist "%TEMP%\gximagecomputing_sunpy" mkdir "%TEMP%\gximagecomputing_sunpy"
if not exist "%TEMP%\gximagecomputing_mpl" mkdir "%TEMP%\gximagecomputing_mpl"

cd /d "%REPO_ROOT%"
set "PYTHONPATH=%REPO_ROOT%\src"
set "SUNPY_CONFIGDIR=%TEMP%\gximagecomputing_sunpy"
set "MPLCONFIGDIR=%TEMP%\gximagecomputing_mpl"

"%PYTHON_CMD%" -c "required = ['gxrender.utils.test_data', 'h5py', 'numpy', 'sunpy.map', 'matplotlib.pyplot']; [__import__(name) for name in required]" >nul 2>&1
if errorlevel 1 (
  echo Failed to import gxrender launcher dependencies with %PYTHON_CMD%.
  exit /b 1
)

for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data root`) do set "TESTDATA_ROOT=%%A"
if not defined TESTDATA_ROOT (
  echo Failed to resolve test-data root.
  exit /b 1
)

if defined MODEL_NAME (
  for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data model "%MODEL_NAME%"`) do set "MODEL=%%A"
) else (
  for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data default-model --suffix .h5`) do set "MODEL=%%A"
)
if not defined MODEL (
  echo Failed to resolve model test data.
  exit /b 1
)

for %%I in ("%MODEL%") do set "MODEL_BASENAME=%%~nxI"
if not defined OUTNAME set "OUTNAME=%MODEL_BASENAME%_py_mw_maps.h5"

for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data ebtel "%EBTEL_NAME%"`) do set "EBTEL=%%A"
if not defined EBTEL (
  echo Failed to resolve EBTEL test data.
  exit /b 1
)

echo Using Python: %PYTHON_CMD%
"%PYTHON_CMD%" src\gxrender\workflows\render_mw.py ^
  --model-path "%MODEL%" ^
  --ebtel-path "%EBTEL%" ^
  --output-dir "%OUTDIR%" ^
  --output-name "%OUTNAME%" ^
  --use-saved-fov
if errorlevel 1 exit /b %errorlevel%

echo You may use gxrender-map-view "%OUTDIR%\%OUTNAME%" to visualize the results
