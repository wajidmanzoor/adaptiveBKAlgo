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

The default binary is configured and incrementally rebuilt in Release mode on
every invocation, so source changes cannot be benchmarked with a stale
executable. Existing output files are never overwritten unless RESUME=1 is
set. In resume mode, completed graph names already present in the CSV are
skipped and new results are appended.

Optional environment variables:
  REORDER_BUDGET       Seed-solver work budget (default: 1000)
  MIN_CLIQUE_SIZE      Minimum reported clique size (default: 3)
  RESUME               Append to OUTPUT_CSV and skip completed rows (0 or 1;
                       default: 0)
  SKIP_GRAPHS          Comma-separated graph basenames to skip; the optional
                       .txt suffix is ignored (default: empty)
  CONTINUE_ON_ERROR    Run later graphs after a graph fails (0 or 1; default: 0)
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
resume=${RESUME:-0}
skip_graphs=${SKIP_GRAPHS:-}
continue_on_error=${CONTINUE_ON_ERROR:-0}
csv_header='graph_name,total_cliques,overall_runtime_seconds'

if [[ ! $reorder_budget =~ ^[0-9]+$ ]]; then
  echo "REORDER_BUDGET must be a nonnegative integer: $reorder_budget" >&2
  exit 2
fi
if [[ ! $minimum_clique_size =~ ^[1-9][0-9]*$ ]]; then
  echo "MIN_CLIQUE_SIZE must be a positive integer: $minimum_clique_size" >&2
  exit 2
fi
if [[ $resume != 0 && $resume != 1 ]]; then
  echo "RESUME must be 0 or 1: $resume" >&2
  exit 2
fi
if [[ $continue_on_error != 0 && $continue_on_error != 1 ]]; then
  echo "CONTINUE_ON_ERROR must be 0 or 1: $continue_on_error" >&2
  exit 2
fi
if [[ ! -d $data_dir ]]; then
  echo "dataset directory does not exist: $data_dir" >&2
  exit 2
fi
if [[ -e $output_csv && $resume == 0 ]]; then
  echo "refusing to overwrite existing output: $output_csv" >&2
  exit 2
fi

if [[ $reorder_binary == "$default_binary" ]]; then
  echo "Configuring and incrementally building optimized Reorder+CCRMCE in Release mode..." >&2
  cmake -S "$reorder_project_dir" -B "$reorder_project_dir/build" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
  cmake --build "$reorder_project_dir/build" --parallel
elif [[ ! -x $reorder_binary ]]; then
  echo "binary is not executable: $reorder_binary" >&2
  exit 2
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

declare -A completed_graphs=()
completed=0
if [[ -e $output_csv ]]; then
  IFS= read -r existing_header < "$output_csv" || true
  if [[ $existing_header != "$csv_header" ]]; then
    echo "cannot resume: unexpected CSV header in $output_csv" >&2
    echo "expected: $csv_header" >&2
    echo "found:    $existing_header" >&2
    exit 2
  fi

  while IFS=, read -r graph_column _; do
    if [[ ${graph_column:0:1} != '"' || ${graph_column: -1} != '"' ]]; then
      echo "cannot resume: malformed graph name in $output_csv: $graph_column" >&2
      exit 2
    fi
    graph_name=${graph_column:1:${#graph_column}-2}
    graph_name=${graph_name//\"\"/\"}
    if [[ -z ${completed_graphs[$graph_name]+present} ]]; then
      completed_graphs["$graph_name"]=1
      (( ++completed ))
    fi
  done < <(sed -n '2,$p' "$output_csv")

  echo "RESUME completed=$completed output=$output_csv" >&2
else
  printf '%s\n' "$csv_header" > "$output_csv"
fi

declare -A requested_skips=()
if [[ -n $skip_graphs ]]; then
  IFS=',' read -r -a skip_entries <<< "$skip_graphs"
  for skip_entry in "${skip_entries[@]}"; do
    if [[ -z $skip_entry ]]; then
      echo "SKIP_GRAPHS contains an empty graph name: $skip_graphs" >&2
      exit 2
    fi
    skip_entry=${skip_entry%.txt}
    requested_skips["$skip_entry"]=1
  done
fi

scratch_dir=$(mktemp -d /tmp/new_reorder_external.XXXXXX)
cleanup() {
  rm -rf -- "$scratch_dir"
}
trap cleanup EXIT

export LC_ALL=C
failures=0
requested_skip_count=0
graph_position=0
total_graphs=${#graph_paths[@]}
for graph_path in "${graph_paths[@]}"; do
  (( ++graph_position ))
  graph_file=$(basename -- "$graph_path")
  graph_name=${graph_file%.txt}
  stdout_file="$scratch_dir/stdout"
  stderr_file="$scratch_dir/stderr"

  if [[ -n ${completed_graphs[$graph_name]+present} ]]; then
    echo "SKIP $graph_position/$total_graphs graph=$graph_file reason=completed" >&2
    continue
  fi
  if [[ -n ${requested_skips[$graph_name]+present} ]]; then
    (( ++requested_skip_count ))
    echo "SKIP $graph_position/$total_graphs graph=$graph_file reason=requested" >&2
    continue
  fi

  echo "RUN $graph_position/$total_graphs graph=$graph_file" >&2
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
    (( ++failures ))
    echo "ERROR graph=$graph_file exit_code=$exit_code" >&2
    sed -n '1,40p' "$stderr_file" >&2
    echo "partial results preserved in: $output_csv" >&2
    if (( continue_on_error == 1 )); then
      echo "CONTINUE graph=$graph_file" >&2
      continue
    fi
    exit "$exit_code"
  fi

  total_cliques=$(awk -F= \
    '$1 == "reorder.cliques" { value=$2 } END { print value }' \
    "$stdout_file")
  if [[ ! $total_cliques =~ ^[0-9]+$ ]]; then
    (( ++failures ))
    echo "ERROR graph=$graph_file did not report reorder.cliques" >&2
    sed -n '1,40p' "$stdout_file" >&2
    echo "partial results preserved in: $output_csv" >&2
    if (( continue_on_error == 1 )); then
      echo "CONTINUE graph=$graph_file" >&2
      continue
    fi
    exit 1
  fi

  escaped_graph_name=${graph_name//\"/\"\"}
  printf '"%s",%s,%s\n' \
    "$escaped_graph_name" "$total_cliques" "$overall_runtime_seconds" \
    >> "$output_csv"
  (( ++completed ))
  echo "DONE graph=$graph_file cliques=$total_cliques runtime_seconds=$overall_runtime_seconds" >&2
done

echo "COMPLETE graphs=$completed/$total_graphs requested_skips=$requested_skip_count failures=$failures output=$output_csv"
if (( failures != 0 )); then
  exit 1
fi
