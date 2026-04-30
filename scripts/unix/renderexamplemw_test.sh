#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTDIR="/tmp/gximagecomputing_validation_groundtruth"
OUTNAME="${OUTNAME:-}"
MODEL_NAME="${MODEL_NAME:-}"
EBTEL_NAME="ebtel.sav"
RUNTIME_CACHE_ROOT="${RUNTIME_CACHE_ROOT:-/tmp/gximagecomputing_runtime_cache}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$RUNTIME_CACHE_ROOT/matplotlib}"
export SUNPY_CONFIGDIR="${SUNPY_CONFIGDIR:-$RUNTIME_CACHE_ROOT/sunpy}"

mkdir -p "$OUTDIR" /tmp/gximagecomputing_sunpy /tmp/gximagecomputing_mpl "$MPLCONFIGDIR" "$SUNPY_CONFIGDIR"

python_supports_render_mw() {
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
    if [[ -x "$CANDIDATE" ]] && python_supports_render_mw "$CANDIDATE"; then
      PYTHON_CMD="$CANDIDATE"
      break
    fi
  done
  if [[ -z "$PYTHON_CMD" ]]; then
    for command_name in python3 python; do
      command_path="$(command -v "$command_name" 2>/dev/null || true)"
      if [[ -n "$command_path" ]] && python_supports_render_mw "$command_path"; then
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

resolve_default_model() {
  resolve_testdata default-model --suffix .h5
}

if [[ -n "$MODEL_NAME" ]]; then
  MODEL="$(resolve_testdata model "$MODEL_NAME")"
else
  MODEL="$(resolve_default_model)"
fi
[[ -n "$MODEL" && -f "$MODEL" ]] || { echo "ERROR: Could not locate an installed H5 model fixture."; exit 1; }
EBTEL="$(resolve_testdata ebtel "$EBTEL_NAME")"
MODEL_BASENAME="$(basename "$MODEL")"
if [[ -z "$OUTNAME" ]]; then
  OUTNAME="${MODEL_BASENAME}_py_mw_maps.h5"
fi

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
# - custom MW frequency list:
#     uncomment `--frequencies-ghz ...`
# - selective heating with default table:
#     uncomment `--selective-heating`
# - selective heating with custom table:
#     uncomment both `--selective-heating` and `--shtable-path ...`
ARGS=(
  --model-path "$MODEL"
  --ebtel-path "$EBTEL"
  --output-dir "$OUTDIR"
  --output-name "$OUTNAME"
  --use-saved-fov

  # Geometry / observer examples
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

  # MW frequency examples
  # --frequencies-ghz 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0 9.4 9.8 10.2 10.6 11.0 11.4 11.8

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

cd "$REPO_ROOT"
echo "Using Python: $PYTHON_CMD"
env PYTHONPATH=src \
    SUNPY_CONFIGDIR="$SUNPY_CONFIGDIR" \
    MPLCONFIGDIR="$MPLCONFIGDIR" \
    "$PYTHON_CMD" src/gxrender/workflows/render_mw.py \
      "${ARGS[@]}"

echo "You may use gxrender-map-view $OUTDIR/$OUTNAME to visualize the results"
