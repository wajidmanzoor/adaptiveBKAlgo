#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_new_reorder_all_external.sh [DATA_DIR [OUTPUT_CSV [BINARY]]]

Runs the optimized Reorder+CCRMCE implementation once on every *.txt graph in
DATA_DIR and writes:

  graph_name,total_cliques,overall_runtime_seconds

Defaults:
  DATA_DIR   /data/labdata/wajid/graphData/adjacencylist/external
  OUTPUT_CSV <project>/results/new_reorder_external.csv
  BINARY     <project>/build/adaptive_bk

The default binary is built in Release mode automatically when it does not
exist. Existing output files are never overwritten.

Optional environment variables:
  REORDER_BUDGET     Seed-solver work budget (default: 1000)
  MIN_CLIQUE_SIZE    Minimum reported clique size (default: 3)
USAGE
}

if (( $# > 3 )); then
  usage >&2
  exit 2
fi
if (( $# == 1 )) && [[ $1 == "-h" || $1 == "--help" ]]; then
  usage
  exit 0
fi

runner_script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
reorder_project_dir=$(cd -- "$runner_script_dir/.." && pwd)

default_data_dir=/data/labdata/wajid/graphData/adjacencylist/external
default_output_csv="$reorder_project_dir/results/new_reorder_external.csv"
default_binary="$reorder_project_dir/build/adaptive_bk"

data_dir=${1:-$default_data_dir}
output_csv=${2:-$default_output_csv}
reorder_binary=${3:-$default_binary}
reorder_budget=${REORDER_BUDGET:-1000}
minimum_clique_size=${MIN_CLIQUE_SIZE:-3}

if [[ ! $reorder_budget =~ ^[0-9]+$ ]]; then
  echo "REORDER_BUDGET must be a nonnegative integer: $reorder_budget" >&2
  exit 2
fi
if [[ ! $minimum_clique_size =~ ^[1-9][0-9]*$ ]]; then
  echo "MIN_CLIQUE_SIZE must be a positive integer: $minimum_clique_size" >&2
  exit 2
fi
if [[ ! -d $data_dir ]]; then
  echo "dataset directory does not exist: $data_dir" >&2
  exit 2
fi
if [[ -e $output_csv ]]; then
  echo "refusing to overwrite existing output: $output_csv" >&2
  exit 2
fi

if [[ ! -x $reorder_binary ]]; then
  if [[ $reorder_binary != "$default_binary" ]]; then
    echo "binary is not executable: $reorder_binary" >&2
    exit 2
  fi
  echo "Building optimized Reorder+CCRMCE in Release mode..." >&2
  cmake -S "$reorder_project_dir" -B "$reorder_project_dir/build" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
  cmake --build "$reorder_project_dir/build" --parallel
fi

mapfile -d '' graph_paths < <(
  find "$data_dir" -maxdepth 1 -type f -name '*.txt' -print0 | sort -z
)
if (( ${#graph_paths[@]} == 0 )); then
  echo "no *.txt graph files found in: $data_dir" >&2
  exit 2
fi

output_dir=$(dirname -- "$output_csv")
mkdir -p -- "$output_dir"
scratch_dir=$(mktemp -d /tmp/new_reorder_external.XXXXXX)
cleanup() {
  rm -rf -- "$scratch_dir"
}
trap cleanup EXIT

export LC_ALL=C
printf '%s\n' 'graph_name,total_cliques,overall_runtime_seconds' > "$output_csv"

completed=0
total_graphs=${#graph_paths[@]}
for graph_path in "${graph_paths[@]}"; do
  graph_file=$(basename -- "$graph_path")
  graph_name=${graph_file%.txt}
  stdout_file="$scratch_dir/stdout"
  stderr_file="$scratch_dir/stderr"

  echo "RUN $((completed + 1))/$total_graphs graph=$graph_file" >&2
  start_nanoseconds=$(date +%s%N)
  "$reorder_binary" "$graph_path" \
    --budget "$reorder_budget" \
    --min-clique-size "$minimum_clique_size" \
    > "$stdout_file" 2> "$stderr_file"
  exit_code=$?
  stop_nanoseconds=$(date +%s%N)
  overall_runtime_seconds=$(awk \
    -v start="$start_nanoseconds" -v stop="$stop_nanoseconds" \
    'BEGIN { printf "%.6f", (stop - start) / 1000000000 }')

  if (( exit_code != 0 )); then
    echo "ERROR graph=$graph_file exit_code=$exit_code" >&2
    sed -n '1,40p' "$stderr_file" >&2
    echo "partial results preserved in: $output_csv" >&2
    exit "$exit_code"
  fi

  total_cliques=$(awk -F= \
    '$1 == "reorder.cliques" { value=$2 } END { print value }' \
    "$stdout_file")
  if [[ ! $total_cliques =~ ^[0-9]+$ ]]; then
    echo "ERROR graph=$graph_file did not report reorder.cliques" >&2
    sed -n '1,40p' "$stdout_file" >&2
    echo "partial results preserved in: $output_csv" >&2
    exit 1
  fi

  escaped_graph_name=${graph_name//\"/\"\"}
  printf '"%s",%s,%s\n' \
    "$escaped_graph_name" "$total_cliques" "$overall_runtime_seconds" \
    >> "$output_csv"
  (( ++completed ))
  echo "DONE graph=$graph_file cliques=$total_cliques runtime_seconds=$overall_runtime_seconds" >&2
done

echo "COMPLETE graphs=$completed output=$output_csv"
