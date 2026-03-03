#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTDIR="/tmp/gximagecomputing_validation_groundtruth"
MODEL="$REPO_ROOT/test_data/test.chr.sav"
RESPONSE="$REPO_ROOT/test_data/aia_response.sav"
EBTEL="/Users/gelu/ssw/packages/gx_simulator/euv/ebtel/ebtel.sav"

mkdir -p "$OUTDIR" /tmp/gximagecomputing_sunpy /tmp/gximagecomputing_mpl

cd "$REPO_ROOT"
env PYTHONPATH=src \
    SUNPY_CONFIGDIR=/tmp/gximagecomputing_sunpy \
    MPLCONFIGDIR=/tmp/gximagecomputing_mpl \
    python src/gximagecomputing/workflows/render_euv.py \
      --model-path "$MODEL" \
      --response-sav "$RESPONSE" \
      --ebtel-path "$EBTEL" \
      --output-dir "$OUTDIR" \
      --output-name test.chr.sav_py_euv_maps.h5

echo "You may use gxrender-map-view /tmp/gximagecomputing_validation_groundtruth/test.chr.sav_py_euv_maps.h5 to visualize the results"
