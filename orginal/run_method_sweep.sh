#!/usr/bin/env bash

set -u
set -o pipefail

# Usage:
#   ./run_method_sweep.sh <graph_dir> [output_csv] [binary] [timeout_seconds]
#
# Example:
#   ./run_method_sweep.sh data sweep.csv build/bk_algorithm 60
#
# The output is a long-format CSV with one row per graph/method:
#   name,category,n,m,mode,method,count,max,checks,time_ms,status
#
# By default this runs the same method set used in the full sweep:
#   0  = Pivot original
#   1  = Pivot ascending degeneracy
#   2  = Pivot descending degeneracy
#   6  = ReorderSib brute force
#   7  = ReorderSib backtracking
#   8  = ReorderSib greedy
#   9  = ReorderSib bitmask
#   10 = ReorderSib inclusion-minimal hitting set
#
# Linux notes:
#   - Uses /usr/bin/env bash, not sh.
#   - Uses GNU timeout when available. If timeout is missing, the script still
#     runs, just without per-run time limits.
#   - Handles spaces in graph paths.
#
# Override modes with:
#   MODES="0 6 7 9 10" ./run_method_sweep.sh data

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRAPH_DIR="${1:-$SCRIPT_DIR/data}"
TIMESTAMP="$(date +"%Y%m%d_%H%M%S")"
OUT_CSV="${2:-$SCRIPT_DIR/sweep_${TIMESTAMP}.csv}"
BINARY="${3:-$SCRIPT_DIR/build/bk_algorithm}"
TIMEOUT_SECONDS="${4:-60}"
MODES="${MODES:-0 1 2 6 7 8 9 10}"

if [ ! -d "$GRAPH_DIR" ]; then
  echo "Error: graph directory not found: $GRAPH_DIR" >&2
  exit 1
fi

if [ ! -x "$BINARY" ]; then
  echo "Error: binary not found or not executable: $BINARY" >&2
  echo "Build first, or pass the binary path as the third argument." >&2
  exit 1
fi

TIMEOUT_CMD=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_CMD="timeout"
elif command -v gtimeout >/dev/null 2>&1; then
  TIMEOUT_CMD="gtimeout"
fi

mode_name() {
  case "$1" in
    0) echo "Pivot-Orig" ;;
    1) echo "Pivot-Asc" ;;
    2) echo "Pivot-Desc" ;;
    6) echo "RSib-BF" ;;
    7) echo "RSib-BT" ;;
    8) echo "RSib-Greedy" ;;
    9) echo "RSib-Bitmask" ;;
    10) echo "RSib-InclMin" ;;
    *) echo "Mode-$1" ;;
  esac
}

infer_category() {
  base="$1"
  parent="$2"
  case "$base" in
    gen_tiny_*) echo "tiny" ;;
    gen_small_*) echo "small" ;;
    gen_medium_*) echo "medium" ;;
    gen_big_*) echo "big" ;;
    data_*) echo "repo_data" ;;
    *) echo "$parent" ;;
  esac
}

csv_escape() {
  value="$1"
  value="${value//\"/\"\"}"
  printf '"%s"' "$value"
}

csv_row() {
  first=1
  for value in "$@"; do
    if [ "$first" -eq 0 ]; then
      printf ','
    fi
    csv_escape "$value"
    first=0
  done
  printf '\n'
}

parse_output() {
  awk '
    /Total Maximal Cliques Found:/ { count=$NF }
    /Maximum Clique Size:/         { maxSize=$NF }
    /Total Vertex-Set Checks:/     { checks=$NF }
    /^Time:/                       { timeMs=$2 }
    /ReorderSib:/ {
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^cliques=/) { split($i, a, "="); count=a[2] }
        if ($i ~ /^maxSize=/) { split($i, a, "="); maxSize=a[2] }
        if ($i ~ /^checks=/)  { split($i, a, "="); checks=a[2] }
        if ($i ~ /^time=/)    { split($i, a, "="); timeMs=a[2] }
      }
    }
    END {
      printf "%s|%s|%s|%s\n", count, maxSize, checks, timeMs
    }
  ' "$1"
}

run_method() {
  graph="$1"
  mode="$2"
  out_file="$3"

  if [ -n "$TIMEOUT_CMD" ]; then
    "$TIMEOUT_CMD" "$TIMEOUT_SECONDS" "$BINARY" "$graph" "$mode" > "$out_file" 2>&1
  else
    "$BINARY" "$graph" "$mode" > "$out_file" 2>&1
  fi
}

graph_list() {
  if sort -z </dev/null >/dev/null 2>&1; then
    find "$GRAPH_DIR" -type f -name "*.txt" -print0 | sort -z
  else
    find "$GRAPH_DIR" -type f -name "*.txt" -print0
  fi
}

mkdir -p "$(dirname "$OUT_CSV")"
csv_row "name" "category" "n" "m" "mode" "method" "count" "max" "checks" "time_ms" "status" > "$OUT_CSV"

graph_count=0
row_count=0

while IFS= read -r -d '' graph; do
  graph_count=$((graph_count + 1))
  name="$(basename "$graph")"
  parent="$(basename "$(dirname "$graph")")"
  category="$(infer_category "$name" "$parent")"

  n=""
  m=""
  if read -r n m _ < "$graph"; then
    :
  fi

  for mode in $MODES; do
    method="$(mode_name "$mode")"
    tmp_out="$(mktemp "${TMPDIR:-/tmp}/bk_sweep.XXXXXX")" || exit 1

    run_method "$graph" "$mode" "$tmp_out"
    rc=$?

    status="ok"
    count=""
    max_size=""
    checks=""
    time_ms=""

    if [ "$rc" -eq 124 ]; then
      status="timeout"
    elif [ "$rc" -ne 0 ]; then
      status="error:$rc"
    else
      parsed="$(parse_output "$tmp_out")"
      IFS='|' read -r count max_size checks time_ms <<< "$parsed"

      if [ -z "$count" ] || [ -z "$max_size" ] || [ -z "$checks" ] || [ -z "$time_ms" ]; then
        status="parse_error"
      fi
    fi

    csv_row "$name" "$category" "$n" "$m" "$mode" "$method" \
      "$count" "$max_size" "$checks" "$time_ms" "$status" >> "$OUT_CSV"

    rm -f "$tmp_out"
    row_count=$((row_count + 1))
  done

  if [ $((graph_count % 10)) -eq 0 ]; then
    echo "Processed $graph_count graphs..."
  fi
done < <(graph_list)

echo "Done."
echo "Graphs processed: $graph_count"
echo "Rows written: $row_count"
echo "CSV: $OUT_CSV"
