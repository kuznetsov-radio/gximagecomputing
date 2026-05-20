#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTDIR="${OUTDIR:-/tmp/gximagecomputing_validation_groundtruth}"
OUTNAME="${OUTNAME:-}"
MODEL_NAME="${MODEL_NAME:-}"
MODEL_PATH="${MODEL_PATH:-}"
EBTEL_NAME="${EBTEL_NAME:-ebtel.sav}"
EBTEL_PATH="${EBTEL_PATH:-}"
INSTRUMENT="${INSTRUMENT:-}"
OBSERVER="${OBSERVER:-}"
RESPONSE="${RESPONSE:-}"
AUTO_FOV="${AUTO_FOV:-0}"
USE_SAVED_FOV="${USE_SAVED_FOV:-0}"
SHOW_MAPS="${SHOW_MAPS:-0}"
XC="${XC:-}"
YC="${YC:-}"
DX="${DX:-}"
DY="${DY:-}"
PIXEL_SCALE_ARCSEC="${PIXEL_SCALE_ARCSEC:-}"
NX="${NX:-}"
NY="${NY:-}"
XRANGE_MIN="${XRANGE_MIN:-}"
XRANGE_MAX="${XRANGE_MAX:-}"
YRANGE_MIN="${YRANGE_MIN:-}"
YRANGE_MAX="${YRANGE_MAX:-}"
DSUN_CM="${DSUN_CM:-}"
LONC_DEG="${LONC_DEG:-}"
B0SUN_DEG="${B0SUN_DEG:-}"
Q0_OVERRIDE="${Q0_OVERRIDE:-}"
HEATING_A="${HEATING_A:-}"
HEATING_B="${HEATING_B:-}"
RUNTIME_CACHE_ROOT="${RUNTIME_CACHE_ROOT:-/tmp/gximagecomputing_runtime_cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$RUNTIME_CACHE_ROOT/matplotlib}"
export SUNPY_CONFIGDIR="${SUNPY_CONFIGDIR:-$RUNTIME_CACHE_ROOT/sunpy}"
CHANNELS=()

usage() {
  cat <<'EOF'
Usage: renderexampleeuv_test.sh [options]

Options:
  --observer NAME           Observer override (earth, stereo-a, stereo-b, solar orbiter)
  --xc ARCSEC              Output map center X (arcsec)
  --yc ARCSEC              Output map center Y (arcsec)
  --dx ARCSEC              Output pixel scale X (arcsec/pixel)
  --dy ARCSEC              Output pixel scale Y (arcsec/pixel)
  --pixel-scale-arcsec V   Shorthand for setting both --dx and --dy
  --nx PIXELS              Output width in pixels
  --ny PIXELS              Output height in pixels
  --xrange XMIN XMAX       Output X range in arcsec
  --yrange YMIN YMAX       Output Y range in arcsec
  --dsun-cm CM             Override observer distance in model before rendering
  --lonc-deg DEG           Override model lonC before rendering
  --b0sun-deg DEG          Override model b0Sun before rendering
  --instrument NAME         EUV instrument override
  --channels CH0 [CH1 ...]  Explicit EUV channel list override
  --response PATH           Explicit response SAV path
  --ebtel PATH              Explicit EBTEL table path (.sav)
  --model-path PATH         Explicit model path
  --model-name NAME         Test-data model fixture name
  --ebtel-name NAME         Test-data EBTEL fixture name (default: ebtel.sav)
  --q0 FLOAT                Override the default closed-field heating normalization
  --a FLOAT                 Override the default closed-field heating exponent a
  --b FLOAT                 Override the default closed-field heating exponent b
  --output-dir PATH         Output directory
  --output-name NAME        Output filename
  --auto-fov                Force inscribing-FOV recomputation
  --use-saved-fov           Force saved-FOV preference
  --show-maps               Launch gxrender-map-view on output when rendering completes
  -h, --help                Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --observer)
      [[ $# -ge 2 ]] || { echo "ERROR: --observer requires a value"; exit 1; }
      OBSERVER="$2"
      shift 2
      ;;
    --xc)
      [[ $# -ge 2 ]] || { echo "ERROR: --xc requires a value"; exit 1; }
      XC="$2"
      shift 2
      ;;
    --yc)
      [[ $# -ge 2 ]] || { echo "ERROR: --yc requires a value"; exit 1; }
      YC="$2"
      shift 2
      ;;
    --dx)
      [[ $# -ge 2 ]] || { echo "ERROR: --dx requires a value"; exit 1; }
      DX="$2"
      shift 2
      ;;
    --dy)
      [[ $# -ge 2 ]] || { echo "ERROR: --dy requires a value"; exit 1; }
      DY="$2"
      shift 2
      ;;
    --pixel-scale-arcsec)
      [[ $# -ge 2 ]] || { echo "ERROR: --pixel-scale-arcsec requires a value"; exit 1; }
      PIXEL_SCALE_ARCSEC="$2"
      shift 2
      ;;
    --nx)
      [[ $# -ge 2 ]] || { echo "ERROR: --nx requires a value"; exit 1; }
      NX="$2"
      shift 2
      ;;
    --ny)
      [[ $# -ge 2 ]] || { echo "ERROR: --ny requires a value"; exit 1; }
      NY="$2"
      shift 2
      ;;
    --xrange)
      [[ $# -ge 3 ]] || { echo "ERROR: --xrange requires XMIN XMAX"; exit 1; }
      XRANGE_MIN="$2"
      XRANGE_MAX="$3"
      shift 3
      ;;
    --yrange)
      [[ $# -ge 3 ]] || { echo "ERROR: --yrange requires YMIN YMAX"; exit 1; }
      YRANGE_MIN="$2"
      YRANGE_MAX="$3"
      shift 3
      ;;
    --dsun-cm)
      [[ $# -ge 2 ]] || { echo "ERROR: --dsun-cm requires a value"; exit 1; }
      DSUN_CM="$2"
      shift 2
      ;;
    --lonc-deg)
      [[ $# -ge 2 ]] || { echo "ERROR: --lonc-deg requires a value"; exit 1; }
      LONC_DEG="$2"
      shift 2
      ;;
    --b0sun-deg)
      [[ $# -ge 2 ]] || { echo "ERROR: --b0sun-deg requires a value"; exit 1; }
      B0SUN_DEG="$2"
      shift 2
      ;;
    --instrument)
      [[ $# -ge 2 ]] || { echo "ERROR: --instrument requires a value"; exit 1; }
      INSTRUMENT="$2"
      shift 2
      ;;
    --channels)
      shift
      [[ $# -ge 1 ]] || { echo "ERROR: --channels requires at least one value"; exit 1; }
      while [[ $# -gt 0 && "$1" != --* ]]; do
        CHANNELS+=("$1")
        shift
      done
      [[ ${#CHANNELS[@]} -gt 0 ]] || { echo "ERROR: --channels requires at least one value"; exit 1; }
      ;;
    --response)
      [[ $# -ge 2 ]] || { echo "ERROR: --response requires a value"; exit 1; }
      RESPONSE="$2"
      shift 2
      ;;
    --ebtel|--ebtel-path)
      [[ $# -ge 2 ]] || { echo "ERROR: --ebtel requires a value"; exit 1; }
      EBTEL_PATH="$2"
      shift 2
      ;;
    --model-path)
      [[ $# -ge 2 ]] || { echo "ERROR: --model-path requires a value"; exit 1; }
      MODEL_PATH="$2"
      shift 2
      ;;
    --model-name)
      [[ $# -ge 2 ]] || { echo "ERROR: --model-name requires a value"; exit 1; }
      MODEL_NAME="$2"
      shift 2
      ;;
    --ebtel-name)
      [[ $# -ge 2 ]] || { echo "ERROR: --ebtel-name requires a value"; exit 1; }
      EBTEL_NAME="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || { echo "ERROR: --output-dir requires a value"; exit 1; }
      OUTDIR="$2"
      shift 2
      ;;
    --output-name)
      [[ $# -ge 2 ]] || { echo "ERROR: --output-name requires a value"; exit 1; }
      OUTNAME="$2"
      shift 2
      ;;
    --q0)
      [[ $# -ge 2 ]] || { echo "ERROR: --q0 requires a value"; exit 1; }
      Q0_OVERRIDE="$2"
      shift 2
      ;;
    --a)
      [[ $# -ge 2 ]] || { echo "ERROR: --a requires a value"; exit 1; }
      HEATING_A="$2"
      shift 2
      ;;
    --b)
      [[ $# -ge 2 ]] || { echo "ERROR: --b requires a value"; exit 1; }
      HEATING_B="$2"
      shift 2
      ;;
    --auto-fov)
      AUTO_FOV="1"
      shift
      ;;
    --use-saved-fov)
      USE_SAVED_FOV="1"
      shift
      ;;
    --show-maps)
      SHOW_MAPS="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

mkdir -p "$OUTDIR" /tmp/gximagecomputing_sunpy /tmp/gximagecomputing_mpl "$MPLCONFIGDIR" "$SUNPY_CONFIGDIR"

python_supports_render_euv() {
  local pycmd="$1"
  (
    cd "$REPO_ROOT"
    env PYTHONPATH=src "$pycmd" -c 'required = ["gxrender.utils.test_data", "h5py", "numpy", "sunpy.map", "matplotlib.pyplot"]; [__import__(name) for name in required]' >/dev/null 2>&1
  )
}

PYTHON_CMD="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_CMD" ]]; then
  CANDIDATES=(
    "$HOME/miniforge3/envs/suncast/bin/python"
    "$HOME/miniforge3/bin/python"
  )
  for CANDIDATE in "${CANDIDATES[@]}"; do
    if [[ -x "$CANDIDATE" ]] && python_supports_render_euv "$CANDIDATE"; then
      PYTHON_CMD="$CANDIDATE"
      break
    fi
  done
  if [[ -z "$PYTHON_CMD" ]]; then
    for command_name in python3 python; do
      command_path="$(command -v "$command_name" 2>/dev/null || true)"
      if [[ -n "$command_path" ]] && python_supports_render_euv "$command_path"; then
        PYTHON_CMD="$command_path"
        break
      fi
    done
  fi
fi
[[ -n "$PYTHON_CMD" ]] || { echo "ERROR: Could not find a Python interpreter with the render-example dependency set."; exit 1; }

resolve_testdata() {
  (
    cd "$REPO_ROOT"
    env PYTHONPATH=src "$PYTHON_CMD" -m gxrender.utils.test_data "$@"
  )
}

resolve_loader_parity_model() {
  (
    cd "$REPO_ROOT"
    env PYTHONPATH=src "$PYTHON_CMD" - <<'PY'
from gxrender.utils.test_data import try_find_model_loader_parity_files

pair = try_find_model_loader_parity_files()
print(pair[1] if pair is not None else "")
PY
  )
}

resolve_default_model() {
  local parity_model
  parity_model="$(resolve_loader_parity_model)"
  if [[ -n "$parity_model" && -f "$parity_model" ]]; then
    printf '%s\n' "$parity_model"
    return 0
  fi
  resolve_testdata default-model --suffix .h5
}

if [[ -n "$MODEL_PATH" ]]; then
  MODEL="$MODEL_PATH"
elif [[ -n "$MODEL_NAME" ]]; then
  MODEL="$(resolve_testdata model "$MODEL_NAME")"
else
  MODEL="$(resolve_default_model)"
fi
[[ -n "$MODEL" && -f "$MODEL" ]] || { echo "ERROR: Could not locate an installed H5 model fixture."; exit 1; }
if [[ -n "$EBTEL_PATH" ]]; then
  EBTEL="$EBTEL_PATH"
else
  EBTEL="$(resolve_testdata ebtel "$EBTEL_NAME")"
fi
[[ -n "$EBTEL" && -f "$EBTEL" ]] || { echo "ERROR: EBTEL file not found: $EBTEL"; exit 1; }
MODEL_BASENAME="$(basename "$MODEL")"
if [[ -z "$OUTNAME" ]]; then
  OUTNAME="${MODEL_BASENAME}_py_euv_maps.h5"
fi
if [[ -n "$RESPONSE" ]]; then
  [[ -f "$RESPONSE" ]] || { echo "ERROR: Response SAV file not found: $RESPONSE"; exit 1; }
fi

MODEL_OBSERVER=""
MODEL_HAS_SAVED_FOV=0
EFFECTIVE_OBSERVER=""
NORMALIZED_OBSERVER=""
NORMALIZED_INSTRUMENT=""
IMPLICIT_RESPONSE_MODE="none"
POLICY_OUTPUT=""
if ! POLICY_OUTPUT="$(
  cd "$REPO_ROOT"
  env PYTHONPATH=src "$PYTHON_CMD" -m gxrender.policy.script_policy \
    --model-path "$MODEL" \
    --instrument "$INSTRUMENT" \
    --observer "$OBSERVER" \
    --response-sav "$RESPONSE"
)"; then
  echo "$POLICY_OUTPUT" >&2
  exit 1
fi

while IFS='=' read -r key value; do
  case "$key" in
    MODEL_OBSERVER) MODEL_OBSERVER="$value" ;;
    HAS_SAVED_FOV) MODEL_HAS_SAVED_FOV="$value" ;;
    EFFECTIVE_OBSERVER) EFFECTIVE_OBSERVER="$value" ;;
    NORMALIZED_OBSERVER) NORMALIZED_OBSERVER="$value" ;;
    NORMALIZED_INSTRUMENT) NORMALIZED_INSTRUMENT="$value" ;;
    IMPLICIT_RESPONSE_MODE) IMPLICIT_RESPONSE_MODE="$value" ;;
  esac
done <<< "$POLICY_OUTPUT"

if [[ "$AUTO_FOV" == "1" && "$USE_SAVED_FOV" == "1" ]]; then
  echo "ERROR: AUTO_FOV=1 and USE_SAVED_FOV=1 are mutually exclusive."
  exit 1
fi

# Top-level fixture selectors:
# - set `MODEL_NAME` to another installed model fixture
# - set `EBTEL_NAME` to another installed EBTEL table
# - set `INSTRUMENT` to force a rendered EUV instrument override
# - set `CHANNELS=(171 193 211)` to force explicit EUV channels
# - set `OBSERVER` to force an observer override (earth, stereo-a, stereo-b, solar orbiter)
# - set `RESPONSE` to force a compatibility SAV override
# - set `AUTO_FOV=1` to force inscribing-FOV recomputation
# - set `USE_SAVED_FOV=1` to force saved-FOV preference
#
# Edit this array to exercise different workflow scenarios.
# Each CLI option is on its own line so you can comment/uncomment it directly.
#
# Suggested patterns:
# - saved observer + saved FOV:
#     keep both `--auto-fov` and `--use-saved-fov` commented out
# - saved observer + recomputed FOV:
#     uncomment `--auto-fov`
# - explicit observer override:
#     uncomment `--observer stereo-a`
# - explicit fixed map box:
#     uncomment `--xc/--yc/--dx/--dy/--nx/--ny`
# - selective heating with default table:
#     uncomment `--selective-heating`
# - selective heating with custom table:
#     uncomment both `--selective-heating` and `--shtable-path ...`
ARGS=(
  --model-path "$MODEL"
  --ebtel-path "$EBTEL"
  --output-dir "$OUTDIR"
  --output-name "$OUTNAME"

  # Geometry / observer examples
  # --use-saved-fov
  # --auto-fov
  # --observer stereo-a
  # --xc -903.0
  # --yc -171.0
  # --dx 2.0
  # --dy 2.0
  # --pixel-scale-arcsec 2.0
  # --nx 75
  # --ny 75
  # --xrange -978.0 -828.0
  # --yrange -246.0 -96.0
  # --dsun-cm 1.4469448e13
  # --lonc-deg -66.84863815759013
  # --b0sun-deg -4.657776399560966

  # Response / channel examples
  # --channels 94 131 171 193 211 304 335

  # Plasma / heating examples
  # --tbase 1.0e6
  # --nbase 1.0e8
  # --q0 0.0217
  # --a 0.3
  # --b 2.7
  # --corona-mode 0
  # --force-isothermal
  # --interpol-b
  # --analytical-nt
  # --selective-heating
  # --shtable-path /path/to/example_shtable.npy
)

if [[ -n "$INSTRUMENT" ]]; then
  ARGS+=(--instrument "$INSTRUMENT")
fi

if [[ -n "$OBSERVER" ]]; then
  ARGS+=(--observer "$OBSERVER")
fi

if [[ -n "$XC" ]]; then
  ARGS+=(--xc "$XC")
fi

if [[ -n "$YC" ]]; then
  ARGS+=(--yc "$YC")
fi

if [[ -n "$DX" ]]; then
  ARGS+=(--dx "$DX")
fi

if [[ -n "$DY" ]]; then
  ARGS+=(--dy "$DY")
fi

if [[ -n "$PIXEL_SCALE_ARCSEC" ]]; then
  ARGS+=(--pixel-scale-arcsec "$PIXEL_SCALE_ARCSEC")
fi

if [[ -n "$NX" ]]; then
  ARGS+=(--nx "$NX")
fi

if [[ -n "$NY" ]]; then
  ARGS+=(--ny "$NY")
fi

if [[ -n "$XRANGE_MIN" && -n "$XRANGE_MAX" ]]; then
  ARGS+=(--xrange "$XRANGE_MIN" "$XRANGE_MAX")
fi

if [[ -n "$YRANGE_MIN" && -n "$YRANGE_MAX" ]]; then
  ARGS+=(--yrange "$YRANGE_MIN" "$YRANGE_MAX")
fi

if [[ -n "$DSUN_CM" ]]; then
  ARGS+=(--dsun-cm "$DSUN_CM")
fi

if [[ -n "$LONC_DEG" ]]; then
  ARGS+=(--lonc-deg "$LONC_DEG")
fi

if [[ -n "$B0SUN_DEG" ]]; then
  ARGS+=(--b0sun-deg "$B0SUN_DEG")
fi

if [[ -n "$RESPONSE" ]]; then
  ARGS+=(--response-sav "$RESPONSE")
fi

if [[ ${#CHANNELS[@]} -gt 0 ]]; then
  ARGS+=(--channels "${CHANNELS[@]}")
fi

if [[ "$AUTO_FOV" == "1" ]]; then
  ARGS+=(--auto-fov)
fi

if [[ "$USE_SAVED_FOV" == "1" ]]; then
  ARGS+=(--use-saved-fov)
fi

if [[ -n "$Q0_OVERRIDE" ]]; then
  ARGS+=(--q0 "$Q0_OVERRIDE")
fi

if [[ -n "$HEATING_A" ]]; then
  ARGS+=(--a "$HEATING_A")
fi

if [[ -n "$HEATING_B" ]]; then
  ARGS+=(--b "$HEATING_B")
fi

cd "$REPO_ROOT"
echo "Using Python: $PYTHON_CMD"
echo "Model observer metadata: ${MODEL_OBSERVER:-<none>}"
echo "Saved FOV metadata present: ${MODEL_HAS_SAVED_FOV}"
if [[ "$IMPLICIT_RESPONSE_MODE" == "earth_aia_time_dependent" ]]; then
  echo "Response selection: implicit Earth/AIA path (time-dependent AIA response inferred from model obs_time)."
fi
env PYTHONPATH=src \
    SUNPY_CONFIGDIR="$SUNPY_CONFIGDIR" \
    MPLCONFIGDIR="$MPLCONFIGDIR" \
    "$PYTHON_CMD" src/gxrender/workflows/render_euv.py \
      "${ARGS[@]}"


echo "You may use gxrender-map-view $OUTDIR/$OUTNAME to visualize the results"
if [[ "$SHOW_MAPS" == "1" ]]; then
  if command -v gxrender-map-view >/dev/null 2>&1; then
    gxrender-map-view "$OUTDIR/$OUTNAME"
  else
    echo "WARNING: gxrender-map-view not found in PATH; cannot auto-open map viewer." >&2
  fi
fi
