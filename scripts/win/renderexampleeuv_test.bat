@echo off
setlocal EnableExtensions EnableDelayedExpansion

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..") do set "REPO_ROOT=%%~fI"

if not defined OUTDIR set "OUTDIR=%TEMP%\gximagecomputing_validation_groundtruth"
if not defined OUTNAME set "OUTNAME="
if not defined MODEL_NAME set "MODEL_NAME="
if not defined MODEL_PATH set "MODEL_PATH="
if not defined EBTEL_NAME set "EBTEL_NAME=ebtel.sav"
if not defined EBTEL_PATH set "EBTEL_PATH="
if not defined INSTRUMENT set "INSTRUMENT="
if not defined OBSERVER set "OBSERVER="
if not defined RESPONSE_SAV if defined RESPONSE set "RESPONSE_SAV=%RESPONSE%"
if not defined RESPONSE_SAV set "RESPONSE_SAV="
if not defined AUTO_FOV set "AUTO_FOV=0"
if not defined USE_SAVED_FOV set "USE_SAVED_FOV=0"
if not defined SHOW_MAPS set "SHOW_MAPS=0"
if not defined XC set "XC="
if not defined YC set "YC="
if not defined DX set "DX="
if not defined DY set "DY="
if not defined PIXEL_SCALE_ARCSEC set "PIXEL_SCALE_ARCSEC="
if not defined NX set "NX="
if not defined NY set "NY="
if not defined XRANGE_MIN set "XRANGE_MIN="
if not defined XRANGE_MAX set "XRANGE_MAX="
if not defined YRANGE_MIN set "YRANGE_MIN="
if not defined YRANGE_MAX set "YRANGE_MAX="
if not defined DSUN_CM set "DSUN_CM="
if not defined LONC_DEG set "LONC_DEG="
if not defined B0SUN_DEG set "B0SUN_DEG="
if not defined CHANNELS set "CHANNELS="
if not defined Q0_OVERRIDE set "Q0_OVERRIDE="
if not defined HEATING_A set "HEATING_A="
if not defined HEATING_B set "HEATING_B="
if not defined RUNTIME_CACHE_ROOT set "RUNTIME_CACHE_ROOT=%TEMP%\gximagecomputing_runtime_cache"
if not defined MPLCONFIGDIR set "MPLCONFIGDIR=%RUNTIME_CACHE_ROOT%\matplotlib"
if not defined SUNPY_CONFIGDIR set "SUNPY_CONFIGDIR=%RUNTIME_CACHE_ROOT%\sunpy"
if defined PYTHON_BIN (
  set "PYTHON_CMD=%PYTHON_BIN%"
) else (
  set "PYTHON_CMD=python"
)

:parse_args
if "%~1"=="" goto args_done
if /i "%~1"=="--observer" (
  if "%~2"=="" (
    echo ERROR: --observer requires a value
    exit /b 1
  )
  set "OBSERVER=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--xc" (
  if "%~2"=="" (
    echo ERROR: --xc requires a value
    exit /b 1
  )
  set "XC=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--yc" (
  if "%~2"=="" (
    echo ERROR: --yc requires a value
    exit /b 1
  )
  set "YC=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--dx" (
  if "%~2"=="" (
    echo ERROR: --dx requires a value
    exit /b 1
  )
  set "DX=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--dy" (
  if "%~2"=="" (
    echo ERROR: --dy requires a value
    exit /b 1
  )
  set "DY=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--pixel-scale-arcsec" (
  if "%~2"=="" (
    echo ERROR: --pixel-scale-arcsec requires a value
    exit /b 1
  )
  set "PIXEL_SCALE_ARCSEC=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--nx" (
  if "%~2"=="" (
    echo ERROR: --nx requires a value
    exit /b 1
  )
  set "NX=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--ny" (
  if "%~2"=="" (
    echo ERROR: --ny requires a value
    exit /b 1
  )
  set "NY=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--xrange" (
  if "%~3"=="" (
    echo ERROR: --xrange requires XMIN XMAX
    exit /b 1
  )
  set "XRANGE_MIN=%~2"
  set "XRANGE_MAX=%~3"
  shift
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--yrange" (
  if "%~3"=="" (
    echo ERROR: --yrange requires YMIN YMAX
    exit /b 1
  )
  set "YRANGE_MIN=%~2"
  set "YRANGE_MAX=%~3"
  shift
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--dsun-cm" (
  if "%~2"=="" (
    echo ERROR: --dsun-cm requires a value
    exit /b 1
  )
  set "DSUN_CM=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--lonc-deg" (
  if "%~2"=="" (
    echo ERROR: --lonc-deg requires a value
    exit /b 1
  )
  set "LONC_DEG=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--b0sun-deg" (
  if "%~2"=="" (
    echo ERROR: --b0sun-deg requires a value
    exit /b 1
  )
  set "B0SUN_DEG=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--instrument" (
  if "%~2"=="" (
    echo ERROR: --instrument requires a value
    exit /b 1
  )
  set "INSTRUMENT=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--channels" goto parse_channels
if /i "%~1"=="--response" (
  if "%~2"=="" (
    echo ERROR: --response requires a value
    exit /b 1
  )
  set "RESPONSE_SAV=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--response-sav" (
  if "%~2"=="" (
    echo ERROR: --response-sav requires a value
    exit /b 1
  )
  set "RESPONSE_SAV=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--ebtel" (
  if "%~2"=="" (
    echo ERROR: --ebtel requires a value
    exit /b 1
  )
  set "EBTEL_PATH=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--ebtel-path" (
  if "%~2"=="" (
    echo ERROR: --ebtel-path requires a value
    exit /b 1
  )
  set "EBTEL_PATH=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--model-path" (
  if "%~2"=="" (
    echo ERROR: --model-path requires a value
    exit /b 1
  )
  set "MODEL_PATH=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--model-name" (
  if "%~2"=="" (
    echo ERROR: --model-name requires a value
    exit /b 1
  )
  set "MODEL_NAME=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--ebtel-name" (
  if "%~2"=="" (
    echo ERROR: --ebtel-name requires a value
    exit /b 1
  )
  set "EBTEL_NAME=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--output-dir" (
  if "%~2"=="" (
    echo ERROR: --output-dir requires a value
    exit /b 1
  )
  set "OUTDIR=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--output-name" (
  if "%~2"=="" (
    echo ERROR: --output-name requires a value
    exit /b 1
  )
  set "OUTNAME=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--q0" (
  if "%~2"=="" (
    echo ERROR: --q0 requires a value
    exit /b 1
  )
  set "Q0_OVERRIDE=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--a" (
  if "%~2"=="" (
    echo ERROR: --a requires a value
    exit /b 1
  )
  set "HEATING_A=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--b" (
  if "%~2"=="" (
    echo ERROR: --b requires a value
    exit /b 1
  )
  set "HEATING_B=%~2"
  shift
  shift
  goto parse_args
)
if /i "%~1"=="--auto-fov" (
  set "AUTO_FOV=1"
  shift
  goto parse_args
)
if /i "%~1"=="--use-saved-fov" (
  set "USE_SAVED_FOV=1"
  shift
  goto parse_args
)
if /i "%~1"=="--show-maps" (
  set "SHOW_MAPS=1"
  shift
  goto parse_args
)
if /i "%~1"=="-h" goto show_help
if /i "%~1"=="--help" goto show_help

echo ERROR: Unknown option: %~1
goto show_help

:parse_channels
shift
if "%~1"=="" (
  echo ERROR: --channels requires at least one value
  exit /b 1
)
:parse_channels_loop
if "%~1"=="" goto parse_channels_done
if /i "%~1"=="--observer" goto parse_channels_done
if /i "%~1"=="--xc" goto parse_channels_done
if /i "%~1"=="--yc" goto parse_channels_done
if /i "%~1"=="--dx" goto parse_channels_done
if /i "%~1"=="--dy" goto parse_channels_done
if /i "%~1"=="--pixel-scale-arcsec" goto parse_channels_done
if /i "%~1"=="--nx" goto parse_channels_done
if /i "%~1"=="--ny" goto parse_channels_done
if /i "%~1"=="--xrange" goto parse_channels_done
if /i "%~1"=="--yrange" goto parse_channels_done
if /i "%~1"=="--dsun-cm" goto parse_channels_done
if /i "%~1"=="--lonc-deg" goto parse_channels_done
if /i "%~1"=="--b0sun-deg" goto parse_channels_done
if /i "%~1"=="--instrument" goto parse_channels_done
if /i "%~1"=="--response" goto parse_channels_done
if /i "%~1"=="--response-sav" goto parse_channels_done
if /i "%~1"=="--ebtel" goto parse_channels_done
if /i "%~1"=="--ebtel-path" goto parse_channels_done
if /i "%~1"=="--model-path" goto parse_channels_done
if /i "%~1"=="--model-name" goto parse_channels_done
if /i "%~1"=="--ebtel-name" goto parse_channels_done
if /i "%~1"=="--output-dir" goto parse_channels_done
if /i "%~1"=="--output-name" goto parse_channels_done
if /i "%~1"=="--q0" goto parse_channels_done
if /i "%~1"=="--a" goto parse_channels_done
if /i "%~1"=="--b" goto parse_channels_done
if /i "%~1"=="--auto-fov" goto parse_channels_done
if /i "%~1"=="--use-saved-fov" goto parse_channels_done
if /i "%~1"=="--show-maps" goto parse_channels_done
if /i "%~1"=="-h" goto parse_channels_done
if /i "%~1"=="--help" goto parse_channels_done
if defined CHANNELS (
  set "CHANNELS=!CHANNELS! %~1"
) else (
  set "CHANNELS=%~1"
)
shift
goto parse_channels_loop

:parse_channels_done
if not defined CHANNELS (
  echo ERROR: --channels requires at least one value
  exit /b 1
)
goto parse_args

:show_help
echo Usage: renderexampleeuv_test.bat [options]
echo.
echo Options:
echo   --observer NAME
echo   --xc ARCSEC
echo   --yc ARCSEC
echo   --dx ARCSEC
echo   --dy ARCSEC
echo   --pixel-scale-arcsec V
echo   --nx PIXELS
echo   --ny PIXELS
echo   --xrange XMIN XMAX
echo   --yrange YMIN YMAX
echo   --dsun-cm CM
echo   --lonc-deg DEG
echo   --b0sun-deg DEG
echo   --instrument NAME
echo   --channels CH0 [CH1 ...]
echo   --response PATH
echo   --ebtel PATH
echo   --model-path PATH
echo   --model-name NAME
echo   --ebtel-name NAME
echo   --q0 FLOAT
echo   --a FLOAT
echo   --b FLOAT
echo   --output-dir PATH
echo   --output-name NAME
echo   --auto-fov
echo   --use-saved-fov
echo   --show-maps
echo   -h, --help
if /i "%~1"=="-h" exit /b 0
if /i "%~1"=="--help" exit /b 0
exit /b 1

:args_done

if not exist "%OUTDIR%" mkdir "%OUTDIR%"
if not exist "%MPLCONFIGDIR%" mkdir "%MPLCONFIGDIR%"
if not exist "%SUNPY_CONFIGDIR%" mkdir "%SUNPY_CONFIGDIR%"

cd /d "%REPO_ROOT%"
set "PYTHONPATH=%REPO_ROOT%\src"

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

if defined MODEL_PATH (
  set "MODEL=%MODEL_PATH%"
) else if defined MODEL_NAME (
  for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data model "%MODEL_NAME%"`) do set "MODEL=%%A"
) else (
  for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -c "from gxrender.utils.test_data import try_find_model_loader_parity_files; pair = try_find_model_loader_parity_files^(^); print^(pair[1] if pair is not None else ''^)"`) do set "MODEL=%%A"
  if not defined MODEL for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data default-model --suffix .h5`) do set "MODEL=%%A"
)
if not defined MODEL (
  echo ERROR: Could not locate an installed H5 model fixture.
  exit /b 1
)
if not exist "%MODEL%" (
  echo ERROR: Could not locate an installed H5 model fixture.
  exit /b 1
)

for %%I in ("%MODEL%") do set "MODEL_BASENAME=%%~nxI"
if not defined OUTNAME set "OUTNAME=%MODEL_BASENAME%_py_euv_maps.h5"
if defined RESPONSE_SAV if not exist "%RESPONSE_SAV%" (
  echo ERROR: RESPONSE_SAV file not found: %RESPONSE_SAV%
  exit /b 1
)

if defined EBTEL_PATH (
  set "EBTEL=%EBTEL_PATH%"
) else (
  for /f "usebackq delims=" %%A in (`"%PYTHON_CMD%" -m gxrender.utils.test_data ebtel "%EBTEL_NAME%"`) do set "EBTEL=%%A"
)
if not defined EBTEL (
  echo Failed to resolve EBTEL test data.
  exit /b 1
)
if not exist "%EBTEL%" (
  echo ERROR: EBTEL file not found: %EBTEL%
  exit /b 1
)

for /f "usebackq tokens=1,* delims==" %%A in (`"%PYTHON_CMD%" -m gxrender.policy.script_policy --model-path "%MODEL%" --instrument "%INSTRUMENT%" --observer "%OBSERVER%" --response-sav "%RESPONSE_SAV%"`) do (
  if /i "%%A"=="MODEL_OBSERVER" set "MODEL_OBSERVER=%%B"
  if /i "%%A"=="HAS_SAVED_FOV" set "MODEL_HAS_SAVED_FOV=%%B"
  if /i "%%A"=="EFFECTIVE_OBSERVER" set "EFFECTIVE_OBSERVER=%%B"
  if /i "%%A"=="NORMALIZED_OBSERVER" set "NORMALIZED_OBSERVER=%%B"
  if /i "%%A"=="NORMALIZED_INSTRUMENT" set "NORMALIZED_INSTRUMENT=%%B"
  if /i "%%A"=="IMPLICIT_RESPONSE_MODE" set "IMPLICIT_RESPONSE_MODE=%%B"
)
if errorlevel 1 exit /b %errorlevel%
if not defined MODEL_HAS_SAVED_FOV set "MODEL_HAS_SAVED_FOV=0"
if not defined IMPLICIT_RESPONSE_MODE set "IMPLICIT_RESPONSE_MODE=none"

if "%AUTO_FOV%"=="1" if "%USE_SAVED_FOV%"=="1" (
  echo ERROR: AUTO_FOV=1 and USE_SAVED_FOV=1 are mutually exclusive.
  exit /b 1
)

set "OBS_ARG="
if defined OBSERVER set "OBS_ARG=--observer ""%OBSERVER%"""

set "GEOM_ARGS="
if defined XC set "GEOM_ARGS=!GEOM_ARGS! --xc !XC!"
if defined YC set "GEOM_ARGS=!GEOM_ARGS! --yc !YC!"
if defined DX set "GEOM_ARGS=!GEOM_ARGS! --dx !DX!"
if defined DY set "GEOM_ARGS=!GEOM_ARGS! --dy !DY!"
if defined PIXEL_SCALE_ARCSEC set "GEOM_ARGS=!GEOM_ARGS! --pixel-scale-arcsec !PIXEL_SCALE_ARCSEC!"
if defined NX set "GEOM_ARGS=!GEOM_ARGS! --nx !NX!"
if defined NY set "GEOM_ARGS=!GEOM_ARGS! --ny !NY!"
if defined XRANGE_MIN if defined XRANGE_MAX set "GEOM_ARGS=!GEOM_ARGS! --xrange !XRANGE_MIN! !XRANGE_MAX!"
if defined YRANGE_MIN if defined YRANGE_MAX set "GEOM_ARGS=!GEOM_ARGS! --yrange !YRANGE_MIN! !YRANGE_MAX!"
if defined DSUN_CM set "GEOM_ARGS=!GEOM_ARGS! --dsun-cm !DSUN_CM!"
if defined LONC_DEG set "GEOM_ARGS=!GEOM_ARGS! --lonc-deg !LONC_DEG!"
if defined B0SUN_DEG set "GEOM_ARGS=!GEOM_ARGS! --b0sun-deg !B0SUN_DEG!"
set "FOV_ARGS="
if "%AUTO_FOV%"=="1" set "FOV_ARGS=!FOV_ARGS! --auto-fov"
if "%USE_SAVED_FOV%"=="1" set "FOV_ARGS=!FOV_ARGS! --use-saved-fov"
set "CHANNEL_ARGS="
if defined CHANNELS set "CHANNEL_ARGS=--channels !CHANNELS!"
set "PLASMA_ARGS="
if defined Q0_OVERRIDE set "PLASMA_ARGS=!PLASMA_ARGS! --q0 !Q0_OVERRIDE!"
if defined HEATING_A set "PLASMA_ARGS=!PLASMA_ARGS! --a !HEATING_A!"
if defined HEATING_B set "PLASMA_ARGS=!PLASMA_ARGS! --b !HEATING_B!"

echo Using Python: %PYTHON_CMD%
echo Model observer metadata: !MODEL_OBSERVER!
echo Saved FOV metadata present: !MODEL_HAS_SAVED_FOV!
if /i "!IMPLICIT_RESPONSE_MODE!"=="earth_aia_time_dependent" (
  echo Response selection: implicit Earth/AIA path (time-dependent AIA response inferred from model obs_time).
)
if defined RESPONSE_SAV (
  "%PYTHON_CMD%" src\gxrender\workflows\render_euv.py ^
    --model-path "%MODEL%" ^
    --response-sav "%RESPONSE_SAV%" ^
    --ebtel-path "%EBTEL%" ^
    --output-dir "%OUTDIR%" ^
    --output-name "%OUTNAME%" ^
    !OBS_ARG! !GEOM_ARGS! !FOV_ARGS! !CHANNEL_ARGS! !PLASMA_ARGS!
) else (
  if defined INSTRUMENT (
    "%PYTHON_CMD%" src\gxrender\workflows\render_euv.py ^
      --model-path "%MODEL%" ^
      --instrument "%INSTRUMENT%" ^
      --ebtel-path "%EBTEL%" ^
      --output-dir "%OUTDIR%" ^
      --output-name "%OUTNAME%" ^
      !OBS_ARG! !GEOM_ARGS! !FOV_ARGS! !CHANNEL_ARGS! !PLASMA_ARGS!
  ) else (
    "%PYTHON_CMD%" src\gxrender\workflows\render_euv.py ^
      --model-path "%MODEL%" ^
      --ebtel-path "%EBTEL%" ^
      --output-dir "%OUTDIR%" ^
      --output-name "%OUTNAME%" ^
      !OBS_ARG! !GEOM_ARGS! !FOV_ARGS! !CHANNEL_ARGS! !PLASMA_ARGS!
  )
)
if errorlevel 1 exit /b %errorlevel%

echo You may use gxrender-map-view "%OUTDIR%\%OUTNAME%" to visualize the results
if "%SHOW_MAPS%"=="1" (
  where gxrender-map-view >nul 2>&1
  if errorlevel 1 (
    echo WARNING: gxrender-map-view not found in PATH; cannot auto-open map viewer. 1>&2
  ) else (
    gxrender-map-view "%OUTDIR%\%OUTNAME%"
  )
)
exit /b 0
