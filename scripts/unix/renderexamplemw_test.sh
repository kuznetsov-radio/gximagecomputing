#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTDIR="/tmp/gximagecomputing_validation_groundtruth"
OUTNAME="test.chr.h5_py_mw_maps.h5"

resolve_testdata() {
  env PYTHONPATH=src python -m gxrender.utils.test_data "$@"
}

MODEL="$(cd "$REPO_ROOT" && resolve_testdata model test.chr.h5)"
EBTEL="$(cd "$REPO_ROOT" && resolve_testdata ebtel ebtel.sav)"

mkdir -p "$OUTDIR" /tmp/gximagecomputing_sunpy /tmp/gximagecomputing_mpl

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

  # Geometry / observer examples
  # --auto-fov
  # --use-saved-fov
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
env PYTHONPATH=src \
    SUNPY_CONFIGDIR=/tmp/gximagecomputing_sunpy \
    MPLCONFIGDIR=/tmp/gximagecomputing_mpl \
    python src/gxrender/workflows/render_mw.py \
      "${ARGS[@]}"

echo "You may use gxrender-map-view $OUTDIR/$OUTNAME to visualize the results"
