#!/usr/bin/env bash

set -uo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
final_root=$(cd "$script_dir/../.." && pwd)
# shellcheck source=../../experiment_common.sh
source "$final_root/experiment_common.sh"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_ROOT [RESULT_ROOT]" >&2
  exit 2
fi

data_root=$1
timestamp=$(date +%Y%m%d_%H%M%S)
result_root=${2:-"$script_dir/results/twelve_$timestamp"}
manifest="$script_dir/../baselines/datasets.tsv"
edge_root="$data_root/edgelist/external"
timeout_seconds=${TOMITA_TIMEOUT_SECONDS:-3600}
build_jobs=${TOMITA_BUILD_JOBS:-4}
dataset_filter=${TOMITA_DATASETS:-}

final_require_positive_integer TOMITA_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer TOMITA_BUILD_JOBS "$build_jobs" || exit 2
final_require_tools make g++ timeout /usr/bin/time awk sha256sum uname mkfifo || exit 2

if [[ ! -f $manifest ]]; then
  echo "Dataset manifest is missing: $manifest" >&2
  exit 2
fi
if [[ ! -d $edge_root ]]; then
  echo "Expected normalized edge lists under: $edge_root" >&2
  exit 2
fi

selected_specs=()
while IFS=$'\t' read -r dataset graph n m reference_cliques; do
  [[ $dataset == dataset ]] && continue
  final_selected_dataset "$dataset" "$dataset_filter" || continue
  if ! final_all_uint "$n" "$m" "$reference_cliques"; then
    echo "Invalid numeric field in manifest row: $dataset" >&2
    exit 2
  fi
  edge_input="$edge_root/$graph"
  if [[ ! -f $edge_input ]]; then
    echo "Missing edge list for $dataset: $edge_input" >&2
    exit 2
  fi
  selected_specs+=(
    "$dataset|$graph|$n|$m|$reference_cliques|$edge_input"
  )
done <"$manifest"

if [[ ${#selected_specs[@]} -eq 0 ]]; then
  echo "No datasets selected." >&2
  exit 2
fi

mkdir -p "$result_root/build_logs" "$result_root/build"

echo "BUILD system=tomita-adjacency-list"
make -C "$script_dir" clean >"$result_root/build_logs/tomita.build.log" 2>&1
make -C "$script_dir" --jobs "$build_jobs" bin/adjlist bin/tomita \
  >>"$result_root/build_logs/tomita.build.log" 2>&1 || {
    echo "Tomita build failed; see $result_root/build_logs/tomita.build.log" >&2
    exit 2
  }
tomita_binary="$script_dir/bin/adjlist"

echo "BUILD helper=tomita_stream_adapter"
adapter="$result_root/build/tomita_stream_adapter"
g++ -O3 -std=c++17 -Wall -Wextra -pedantic \
  "$script_dir/../baselines/rmce_stream_adapter.cpp" -o "$adapter" \
  >"$result_root/build_logs/tomita_stream_adapter.build.log" 2>&1 || {
    echo "Input-adapter build failed; see the build log." >&2
    exit 2
  }

if [[ ! -x $tomita_binary || ! -x $adapter ]]; then
  echo "A benchmark executable is missing after the build." >&2
  exit 2
fi

{
  echo "campaign=Quick Cliques v1.0 Tomita adjacency-list twelve-dataset run"
  echo "upstream_tag=v1.0"
  echo "upstream_commit=dd9c4d9fdf8243e2861aa6b99e8992836da741a9"
  echo "data_root=$data_root"
  echo "edge_root=$edge_root"
  echo "manifest=$manifest"
  echo "result_root=$result_root"
  echo "datasets=${#selected_specs[@]}"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "minimum_clique_size=3"
  echo "clique_storage=all"
  echo "clique_output=disabled"
  echo "algorithm_timer=CLOCK_MONOTONIC; includes explicit clique insertion"
  echo "process_wall_timer=includes adapter, input, enumeration, and cleanup"
  echo "implementation=tomita-adjacency-list (upstream sparse Tomita variant)"
  echo "matrix_tomita_not_used=n-squared storage is infeasible for the full corpus"
  sha256sum "$tomita_binary" "$adapter" "$manifest"
  uname -a
} >"$result_root/environment.txt"

csv="$result_root/results.csv"
if [[ ! -f $csv ]]; then
  printf '%s\n' \
    'system,dataset,graph,n,m,status,exit_code,cliques,stored_cliques,reference_cliques,count_match,algorithm_wall_ms,process_wall_ms,max_rss_kb,input_adapter,adapter_exit_code' \
    >"$csv"
fi

already_recorded() {
  local dataset=$1
  awk -F, -v dataset="$dataset" \
    'NR > 1 && $2 == dataset { found = 1 } END { exit !found }' "$csv"
}

campaign_failed=0

run_dataset() {
  local dataset=$1 graph=$2 n=$3 m=$4 reference_cliques=$5 edge_input=$6
  if already_recorded "$dataset"; then
    echo "SKIP system=tomita-adjacency-list dataset=$dataset reason=recorded"
    return
  fi

  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local log_dir="$result_root/logs/$safe_graph"
  local base="$log_dir/tomita_adjacency_list"
  local fifo="$log_dir/input.fifo"
  mkdir -p "$log_dir"
  if [[ -e $fifo || -p $fifo ]]; then
    rm -f -- "$fifo"
  fi
  mkfifo "$fifo" || {
    echo "Could not create input FIFO: $fifo" >&2
    campaign_failed=1
    return
  }

  final_write_command "$base.command" env OMP_NUM_THREADS=1 \
    "$tomita_binary" "$fifo"
  final_write_command "$base.adapter.command" "$adapter" \
    "$edge_input" "$n" "$m"

  "$adapter" "$edge_input" "$n" "$m" >"$fifo" \
    2>"$base.adapter.stderr" &
  local adapter_pid=$!
  final_run_timed \
    "system=tomita-adjacency-list dataset=$dataset" \
    "$base.stdout" "$base.stderr" "$base.resources" "$timeout_seconds" \
    env OMP_NUM_THREADS=1 "$tomita_binary" "$fifo"

  local exit_code=$FINAL_RUN_EXIT_CODE
  local process_wall_ms=$FINAL_RUN_WALL_MS
  if kill -0 "$adapter_pid" 2>/dev/null; then
    kill "$adapter_pid" 2>/dev/null || true
  fi
  wait "$adapter_pid" 2>/dev/null
  local adapter_exit=$?
  rm -f -- "$fifo"

  local algorithm minimum storage output cliques stored algorithm_wall_ms
  local max_rss_kb status count_match=NA
  algorithm=$(final_output_value "$base.stdout" algorithm)
  minimum=$(final_output_value "$base.stdout" minimum_clique_size)
  storage=$(final_output_value "$base.stdout" clique_storage)
  output=$(final_output_value "$base.stdout" clique_output)
  cliques=$(final_output_value "$base.stdout" maximal_cliques)
  stored=$(final_output_value "$base.stdout" stored_cliques)
  algorithm_wall_ms=$(final_output_value "$base.stdout" algorithm_wall_ms)
  max_rss_kb=$(final_output_value "$base.resources" max_rss_kb)

  if [[ $exit_code -eq 0 && $adapter_exit -eq 0 ]] &&
      final_all_uint "$cliques" "$stored" "$process_wall_ms" "$max_rss_kb" &&
      final_all_number "$algorithm_wall_ms"; then
    if [[ $algorithm != tomita-adjacency-list || $minimum != 3 ||
          $storage != all || $output != disabled ]]; then
      status=config_mismatch
      campaign_failed=1
    elif [[ $cliques != "$stored" ]]; then
      status=storage_mismatch
      campaign_failed=1
    elif [[ $cliques == "$reference_cliques" ]]; then
      status=completed
      count_match=yes
    else
      status=count_mismatch
      count_match=no
      campaign_failed=1
    fi
  else
    status=$(final_status_from_exit "$exit_code")
    campaign_failed=1
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    tomita-adjacency-list "$dataset" "$graph" "$n" "$m" "$status" \
    "$exit_code" "${cliques:-NA}" "${stored:-NA}" "$reference_cliques" \
    "$count_match" "${algorithm_wall_ms:-NA}" "$process_wall_ms" \
    "${max_rss_kb:-NA}" streamed_symmetric_comma_fifo "$adapter_exit" \
    >>"$csv"
  echo "DONE system=tomita-adjacency-list dataset=$dataset status=$status cliques=${cliques:-NA} algorithm_wall_ms=${algorithm_wall_ms:-NA} process_wall_ms=$process_wall_ms"
}

for spec in "${selected_specs[@]}"; do
  IFS='|' read -r dataset graph n m reference_cliques edge_input <<<"$spec"
  run_dataset "$dataset" "$graph" "$n" "$m" "$reference_cliques" \
    "$edge_input"
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${TOMITA_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
