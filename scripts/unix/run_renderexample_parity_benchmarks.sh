#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MODE="${MODE:-both}"
DEFAULT_TMP_ROOT="${RUNNER_TEMP:-${TMPDIR:-/tmp}}"
OUTROOT="${OUTROOT:-$DEFAULT_TMP_ROOT/gximagecomputing_renderexample_parity}"
RUNTIME_CACHE_ROOT="${RUNTIME_CACHE_ROOT:-$DEFAULT_TMP_ROOT/gximagecomputing_runtime_cache}"
PYTHON_CMD="${PYTHON_BIN:-}"
IDL_CMD="${IDL_BIN:-}"

export MPLCONFIGDIR="${MPLCONFIGDIR:-$RUNTIME_CACHE_ROOT/matplotlib}"
export SUNPY_CONFIGDIR="${SUNPY_CONFIGDIR:-$RUNTIME_CACHE_ROOT/sunpy}"

usage() {
  cat <<'EOF'
Usage: run_renderexample_parity_benchmarks.sh [options]

Runs Python-vs-IDL renderexample parity benchmarks for MW, EUV, or both.

Options:
  --mode mw|euv|both       Benchmark mode (default: MODE or both)
  --outroot PATH           Output artifact root (default: OUTROOT or CI temp)
  --python PATH            Python interpreter (default: PYTHON_BIN or PATH lookup)
  --idl PATH               IDL launcher (default: IDL_BIN or PATH lookup)
  --model-h5 PATH          Python workflow H5 model input
  --model-sav PATH         IDL workflow SAV model input
  --model-path PATH        Backward-compatible fallback for both model inputs
  --ebtel PATH             EBTEL SAV fixture path
  --response-sav PATH      EUV response SAV fixture path
  -h, --help               Show this help

Environment equivalents:
  MODE, OUTROOT, PYTHON_BIN, IDL_BIN, MODEL_H5_PATH, MODEL_SAV_PATH,
  MODEL_PATH, GXIMAGECOMPUTING_EBTEL_PATH, GXIMAGECOMPUTING_EUV_RESPONSE_SAV,
  RUNTIME_CACHE_ROOT, MPLCONFIGDIR, SUNPY_CONFIGDIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      [[ $# -ge 2 ]] || { echo "ERROR: --mode requires a value" >&2; exit 1; }
      MODE="$2"
      shift 2
      ;;
    --outroot)
      [[ $# -ge 2 ]] || { echo "ERROR: --outroot requires a value" >&2; exit 1; }
      OUTROOT="$2"
      shift 2
      ;;
    --python)
      [[ $# -ge 2 ]] || { echo "ERROR: --python requires a value" >&2; exit 1; }
      PYTHON_CMD="$2"
      shift 2
      ;;
    --idl)
      [[ $# -ge 2 ]] || { echo "ERROR: --idl requires a value" >&2; exit 1; }
      IDL_CMD="$2"
      shift 2
      ;;
    --model-h5)
      [[ $# -ge 2 ]] || { echo "ERROR: --model-h5 requires a value" >&2; exit 1; }
      MODEL_H5_PATH="$2"
      shift 2
      ;;
    --model-sav)
      [[ $# -ge 2 ]] || { echo "ERROR: --model-sav requires a value" >&2; exit 1; }
      MODEL_SAV_PATH="$2"
      shift 2
      ;;
    --model-path)
      [[ $# -ge 2 ]] || { echo "ERROR: --model-path requires a value" >&2; exit 1; }
      MODEL_PATH="$2"
      shift 2
      ;;
    --ebtel)
      [[ $# -ge 2 ]] || { echo "ERROR: --ebtel requires a value" >&2; exit 1; }
      GXIMAGECOMPUTING_EBTEL_PATH="$2"
      shift 2
      ;;
    --response-sav)
      [[ $# -ge 2 ]] || { echo "ERROR: --response-sav requires a value" >&2; exit 1; }
      GXIMAGECOMPUTING_EUV_RESPONSE_SAV="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

mkdir -p "$OUTROOT" "$MPLCONFIGDIR" "$SUNPY_CONFIGDIR"

resolve_command() {
  local configured="$1"
  shift
  local candidate

  if [[ -n "$configured" ]]; then
    if [[ -x "$configured" ]]; then
      printf '%s\n' "$configured"
      return 0
    fi
    if command -v "$configured" >/dev/null 2>&1; then
      command -v "$configured"
      return 0
    fi
    echo "ERROR: Command not found or not executable: $configured" >&2
    return 1
  fi

  for candidate in "$@"; do
    if command -v "$candidate" >/dev/null 2>&1; then
      command -v "$candidate"
      return 0
    fi
  done
  echo "ERROR: Could not find any of: $*" >&2
  return 1
}

PYTHON_CMD="$(resolve_command "$PYTHON_CMD" python3 python)"
IDL_CMD="$(resolve_command "$IDL_CMD" sswidl idl)"

run_py() {
  (
    cd "$REPO_ROOT"
    env PYTHONPATH=src "$PYTHON_CMD" "$@"
  )
}

resolve_loader_parity_pair() {
  run_py - <<'PY'
from gxrender.utils.test_data import find_model_loader_parity_files
sav_path, h5_path = find_model_loader_parity_files()
print(sav_path)
print(h5_path)
PY
}

resolve_ebtel() {
  run_py -m gxrender.utils.test_data ebtel ebtel.sav
}

resolve_response() {
  run_py -m gxrender.utils.test_data response aia
}

summarize_compare_json() {
  local json_path="$1"
  run_py - "$json_path" <<'PY'
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
mode = report.get('mode', {})
summary = report.get('summary', {})
comp1 = mode.get('comp1_label', 'COMP1')
comp2 = mode.get('comp2_label', 'COMP2')

def show(label: str) -> None:
    abs_stats = summary.get(f"{label}_abs_diff", {})
    pct_stats = summary.get(f"{label}_pct_diff_vs_idl", {})
    sym_stats = summary.get(f"{label}_sym_diff_pm1", {})
    print(f"- {label}: mean_abs_diff={abs_stats.get('mean')}, max_abs_diff={abs_stats.get('max')}, "
          f"mean_pct_diff={pct_stats.get('mean')}, max_pct_diff={pct_stats.get('max')}, "
          f"max_sym_diff={sym_stats.get('max')}")

print(f"- comparison_json: {sys.argv[1]}")
show(comp1)
show(comp2)
PY
}

run_idl_batch() {
  local batch_file="$1"
  "$IDL_CMD" < "$batch_file"
}

make_temp_batch() {
  mktemp "${TMPDIR:-/tmp}/gximage_renderexample_parity.XXXXXX"
}

write_idl_batch() {
  local batch_file="$1"
  local proc_name="$2"
  local model_path="$3"
  local ebtel_path="$4"
  local out_dir="$5"
  local out_file="$6"
  local response_sav="${7:-}"

  {
    printf '%s\n' "@${REPO_ROOT}/examples/idl/compile_local_idl"
    printf '%s' "${proc_name}, MODelfile='${model_path}', EBTELfile='${ebtel_path}', OUTdir='${out_dir}', OUTfile='${out_file}', USE_SAVED_FOV=1, NO_PLOT=1"
    if [[ -n "$response_sav" ]]; then
      printf '%s' ", RESPonsefile='${response_sav}'"
    fi
    printf '\nexit\n'
  } > "$batch_file"
}

run_mw_benchmark() {
  local model_h5="$1"
  local model_sav="$2"
  local mw_dir="$OUTROOT/mw"
  local ebtel_path="$3"
  local py_h5="$mw_dir/$(basename "$model_h5")_py_mw_maps.h5"
  local idl_sav="$mw_dir/$(basename "$model_sav")_idl_mw_maps.sav"
  local compare_json="$mw_dir/comparison_python_vs_idl.json"
  local batch_file

  mkdir -p "$mw_dir"
  echo "Running MW Python renderexample benchmark..."
  (
    cd "$REPO_ROOT"
    OUTDIR="$mw_dir" \
    OUTNAME="$(basename "$py_h5")" \
    MODEL_PATH="$model_h5" \
    PYTHON_BIN="$PYTHON_CMD" \
    RUNTIME_CACHE_ROOT="$RUNTIME_CACHE_ROOT" \
    bash scripts/unix/renderexamplemw_test.sh
  )

  echo "Running MW IDL renderexample benchmark..."
  batch_file="$(make_temp_batch)"
  trap 'rm -f "$batch_file"' RETURN
  write_idl_batch "$batch_file" "RenderExampleMW_test" "$model_sav" "$ebtel_path" "$mw_dir" "$(basename "$idl_sav")"
  run_idl_batch "$batch_file"
  rm -f "$batch_file"
  trap - RETURN

  echo "Comparing MW outputs..."
  run_py scripts/python/ComparePythonVsIDLMaps.py --kind mw --out-dir "$mw_dir" --python-h5 "$py_h5" --idl-sav "$idl_sav"
  summarize_compare_json "$compare_json"
}

run_euv_benchmark() {
  local model_h5="$1"
  local model_sav="$2"
  local euv_dir="$OUTROOT/euv"
  local ebtel_path="$3"
  local response_sav="$4"
  local py_h5="$euv_dir/$(basename "$model_h5")_py_euv_maps.h5"
  local idl_sav="$euv_dir/$(basename "$model_sav")_idl_euv_maps.sav"
  local compare_json="$euv_dir/comparison_python_vs_idl.json"
  local batch_file

  mkdir -p "$euv_dir"
  echo "Running EUV Python renderexample benchmark..."
  (
    cd "$REPO_ROOT"
    OUTDIR="$euv_dir" \
    OUTNAME="$(basename "$py_h5")" \
    MODEL_PATH="$model_h5" \
    RESPONSE_SAV="$response_sav" \
    PYTHON_BIN="$PYTHON_CMD" \
    RUNTIME_CACHE_ROOT="$RUNTIME_CACHE_ROOT" \
    bash scripts/unix/renderexampleeuv_test.sh
  )

  echo "Running EUV IDL renderexample benchmark..."
  batch_file="$(make_temp_batch)"
  trap 'rm -f "$batch_file"' RETURN
  write_idl_batch "$batch_file" "RenderExampleEUV_test" "$model_sav" "$ebtel_path" "$euv_dir" "$(basename "$idl_sav")" "$response_sav"
  run_idl_batch "$batch_file"
  rm -f "$batch_file"
  trap - RETURN

  echo "Comparing EUV outputs..."
  run_py scripts/python/ComparePythonVsIDLEUVMaps.py --out-dir "$euv_dir" --python-h5 "$py_h5" --idl-sav "$idl_sav"
  summarize_compare_json "$compare_json"
}

main() {
  local model_sav
  local model_h5
  local ebtel_path
  local response_sav=""
  local pair_output

  case "$MODE" in
    mw|euv|both) ;;
    *)
      echo "ERROR: MODE must be one of mw, euv, both." >&2
      exit 1
      ;;
  esac

  if [[ -n "${MODEL_H5_PATH:-}" || -n "${MODEL_SAV_PATH:-}" ]]; then
    [[ -n "${MODEL_H5_PATH:-}" ]] || { echo "ERROR: MODEL_H5_PATH/--model-h5 is required when MODEL_SAV_PATH is set." >&2; exit 1; }
    [[ -n "${MODEL_SAV_PATH:-}" ]] || { echo "ERROR: MODEL_SAV_PATH/--model-sav is required when MODEL_H5_PATH is set." >&2; exit 1; }
    model_h5="$MODEL_H5_PATH"
    model_sav="$MODEL_SAV_PATH"
  elif [[ -n "${MODEL_PATH:-}" ]]; then
    echo "WARNING: MODEL_PATH/--model-path is using the same file for Python and IDL inputs; prefer MODEL_H5_PATH and MODEL_SAV_PATH in CI." >&2
    model_h5="$MODEL_PATH"
    model_sav="$MODEL_PATH"
  else
    pair_output="$(resolve_loader_parity_pair)"
    model_sav="$(printf '%s\n' "$pair_output" | sed -n '1p')"
    model_h5="$(printf '%s\n' "$pair_output" | sed -n '2p')"
  fi

  [[ -f "$model_h5" ]] || { echo "ERROR: H5 model not found: $model_h5" >&2; exit 1; }
  [[ -f "$model_sav" ]] || { echo "ERROR: SAV model not found: $model_sav" >&2; exit 1; }

  ebtel_path="${GXIMAGECOMPUTING_EBTEL_PATH:-$(resolve_ebtel)}"
  [[ -f "$ebtel_path" ]] || { echo "ERROR: EBTEL file not found: $ebtel_path" >&2; exit 1; }
  if [[ "$MODE" == "euv" || "$MODE" == "both" ]]; then
    response_sav="${GXIMAGECOMPUTING_EUV_RESPONSE_SAV:-$(resolve_response)}"
    [[ -f "$response_sav" ]] || { echo "ERROR: EUV response SAV not found: $response_sav" >&2; exit 1; }
  fi

  echo "Benchmark fixtures:"
  echo "- python: $PYTHON_CMD"
  echo "- idl: $IDL_CMD"
  echo "- model_h5: $model_h5"
  echo "- model_sav: $model_sav"
  echo "- ebtel: $ebtel_path"
  if [[ -n "$response_sav" ]]; then
    echo "- response_sav: $response_sav"
  fi
  echo "- outroot: $OUTROOT"

  case "$MODE" in
    mw)
      run_mw_benchmark "$model_h5" "$model_sav" "$ebtel_path"
      ;;
    euv)
      run_euv_benchmark "$model_h5" "$model_sav" "$ebtel_path" "$response_sav"
      ;;
    both)
      run_mw_benchmark "$model_h5" "$model_sav" "$ebtel_path"
      run_euv_benchmark "$model_h5" "$model_sav" "$ebtel_path" "$response_sav"
      ;;
  esac

  echo "Benchmark artifacts written under: $OUTROOT"
}

main "$@"
