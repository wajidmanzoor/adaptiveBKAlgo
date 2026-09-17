#!/usr/bin/env bash

set -uo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=../../experiment_common.sh
source "$script_dir/../../experiment_common.sh"
# shellcheck source=../../grouped_experiment_common.sh
source "$script_dir/../../grouped_experiment_common.sh"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 DATA_ROOT [RESULT_ROOT]" >&2
  exit 2
fi

data_root=$1
timestamp=$(date +%Y%m%d_%H%M%S)
result_root=${2:-"$script_dir/results/$timestamp"}
timeout_seconds=${CAPACITY_TIMEOUT_SECONDS:-3600}
build_jobs=${CAPACITY_BUILD_JOBS:-4}
dataset_filter=${CAPACITY_DATASETS:-}
repetitions=${CAPACITY_REPETITIONS:-3}
budget=1000

capacities=(64 128 512 1024 dynamic)
# Run the reference first so every subsequent row can be checked immediately.
run_capacities=(128 64 512 1024 dynamic)

final_require_positive_integer CAPACITY_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer CAPACITY_BUILD_JOBS "$build_jobs" || exit 2
final_require_positive_integer CAPACITY_REPETITIONS "$repetitions" || exit 2
final_require_tools cmake timeout /usr/bin/time awk sha256sum uname find sort || exit 2
final_collect_grouped_datasets "$data_root" "$dataset_filter" 0 || exit 2

mkdir -p "$result_root/build_logs"
pure_source="$script_dir/pure"
declare -A binaries
for capacity in "${capacities[@]}"; do
  build_root="$result_root/build/capacity_$capacity"
  final_build "capacity_$capacity" "$pure_source" "$build_root" \
    "$result_root/build_logs" "$build_jobs" \
    -DPURE_HITSET_VARIANT="$capacity" || exit 2
  binaries[$capacity]="$build_root/$FINAL_PURE_BINARY_NAME"
  [[ -x ${binaries[$capacity]} ]] || {
    echo "Missing reorder binary for capacity $capacity." >&2
    exit 2
  }
done

{
  echo "campaign=final Reorder+CCRMCE seed-mask capacity ablation"
  echo "data_root=$data_root"
  echo "adjacency_root=$FINAL_ADJACENCY_ROOT"
  echo "result_root=$result_root"
  echo "groups=${FINAL_SELECTED_GROUPS[*]}"
  echo "graphs=${#FINAL_SELECTED_DATASETS[@]}"
  echo "capacities=${capacities[*]}"
  echo "reference_capacity=128"
  echo "budget=$budget (fixed)"
  echo "minimum_clique_size=3"
  echo "small_q_ccr_threshold=32"
  echo "early_termination=ET1 ET2 ET3 disabled"
  echo "pruning=production profile (subsumption disabled; six other rules enabled)"
  echo "repetitions=$repetitions"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "execution=sequential"
  echo "dynamic_semantics=ceil(reduced_constraints/64) words per nontrivial seed-solver call"
  for capacity in "${capacities[@]}"; do
    sha256sum "${binaries[$capacity]}"
  done
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local runs_csv=$1
  local capacity=$2
  local graph=$3
  local repetition=$4
  awk -F, -v capacity="$capacity" -v graph="$graph" \
    -v repetition="$repetition" \
    'NR > 1 && $1 == capacity && $2 == graph && $3 == repetition {
       found = 1
     }
     END { exit !found }' "$runs_csv"
}

reference_field() {
  local runs_csv=$1
  local graph=$2
  local repetition=$3
  local field=$4
  awk -F, -v graph="$graph" -v repetition="$repetition" -v field="$field" \
    'NR > 1 && $1 == "128" && $2 == graph && $3 == repetition {
       value = $field
     }
     END { print value }' "$runs_csv"
}

run_one() {
  local capacity=$1
  local graph=$2
  local repetition=$3
  local n=$4
  local m=$5
  local input=$6
  local group_root=$7
  local runs_csv=$8

  if already_recorded "$runs_csv" "$capacity" "$graph" "$repetition"; then
    echo "SKIP capacity=$capacity graph=$graph repetition=$repetition reason=recorded"
    return
  fi

  local binary=${binaries[$capacity]}
  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local log_dir="$group_root/logs/$safe_graph"
  local base="$log_dir/capacity_${capacity}__r${repetition}"
  mkdir -p "$log_dir"
  local -a command=(
    env OMP_NUM_THREADS=1 "$binary" "$input"
    --budget "$budget" --min-clique-size 3
  )

  final_write_command "$base.command" "${command[@]}"
  final_run_timed \
    "capacity=$capacity graph=$graph repetition=$repetition" \
    "$base.stdout" "$base.stderr" "$base.resources" "$timeout_seconds" \
    "${command[@]}"

  local exit_code=$FINAL_RUN_EXIT_CODE
  local wall_ms=$FINAL_RUN_WALL_MS
  local cliques stored runtime_ms max_rss findone_states full_states ccr_states budget_fallbacks
  local capacity_fallbacks solver_calls maximum_constraints configured dynamic
  local expected_dynamic status reference_status reference_count reference_wall
  local count_match=NA wall_ratio=NA

  cliques=$(final_output_value "$base.stdout" reorder.cliques)
  stored=$(final_output_value "$base.stdout" reorder.stored_cliques)
  runtime_ms=$(final_output_value "$base.stdout" reorder.runtime_ms)
  max_rss=$(final_output_value "$base.resources" max_rss_kb)
  findone_states=$(final_output_value "$base.stdout" reorder.ccr.findone_states)
  full_states=$(final_output_value "$base.stdout" reorder.ccr.full_states)
  ccr_states=
  if final_all_uint "$findone_states" "$full_states"; then
    ccr_states=$((findone_states + full_states))
  fi
  budget_fallbacks=$(final_output_value "$base.stdout" reorder.budget_fallbacks)
  capacity_fallbacks=$(final_output_value "$base.stdout" reorder.capacity_fallbacks)
  solver_calls=$(final_output_value "$base.stdout" reorder.seed_solver_calls)
  maximum_constraints=$(final_output_value "$base.stdout" reorder.maximum_seed_constraints)
  configured=$(final_output_value "$base.stdout" reorder.config.hitset_capacity)
  dynamic=$(final_output_value "$base.stdout" reorder.config.hitset_dynamic)

  expected_dynamic=0
  [[ $capacity == dynamic ]] && expected_dynamic=1

  if [[ $exit_code -eq 0 ]] &&
      final_all_uint "$cliques" "$stored" "$max_rss" "$ccr_states" \
        "$budget_fallbacks" "$capacity_fallbacks" "$solver_calls" \
        "$maximum_constraints" &&
      final_all_number "$runtime_ms" &&
      [[ $stored == "$cliques" && $configured == "$capacity" &&
         $dynamic == "$expected_dynamic" ]] &&
      [[ $(final_output_value "$base.stdout" reorder.budget) == "$budget" ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.small_q_ccr_threshold) == 32 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et1) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et2) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et3) == 0 ]] &&
      final_pruning_config_matches "$base.stdout" production; then
    status=completed
  else
    status=$(final_status_from_exit "$exit_code")
  fi

  if [[ $capacity == 128 ]]; then
    reference_status=$status
    reference_count=$cliques
    reference_wall=$wall_ms
  else
    reference_status=$(reference_field "$runs_csv" "$graph" "$repetition" 6)
    reference_count=$(reference_field "$runs_csv" "$graph" "$repetition" 8)
    reference_wall=$(reference_field "$runs_csv" "$graph" "$repetition" 11)
  fi

  if [[ $status == completed && $reference_status == completed ]]; then
    [[ $cliques == "$reference_count" ]] && count_match=yes || count_match=no
    if [[ $reference_wall =~ ^[0-9]+$ && $reference_wall -gt 0 ]]; then
      wall_ratio=$(awk -v wall="$wall_ms" -v reference="$reference_wall" \
        'BEGIN { printf "%.6f", wall / reference }')
    fi
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$capacity" "$graph" "$repetition" "$n" "$m" "$status" \
    "$exit_code" "$cliques" "$stored" "$count_match" "$wall_ms" \
    "$wall_ratio" "$runtime_ms" "$max_rss" "$ccr_states" \
    "$budget_fallbacks" "$capacity_fallbacks" "$solver_calls" \
    "$maximum_constraints" "$configured" "$dynamic" "$budget" \
    >>"$runs_csv"

  echo "DONE capacity=$capacity graph=$graph repetition=$repetition status=$status cliques=${cliques:-NA} match=$count_match wall_ratio=$wall_ratio capacity_fallbacks=${capacity_fallbacks:-NA}"
}

write_summary() {
  local runs_csv=$1
  local summary_csv=$2
  printf '%s\n' \
    'capacity,completed_runs,total_runs,paired_runs,geomean_wall_ratio_vs_128,geomean_speedup_vs_128,mean_wall_ms,mean_runtime_ms,mean_max_rss_kb,total_budget_fallbacks,total_capacity_fallbacks,maximum_seed_constraints,all_counts_match' \
    >"$summary_csv"

  local capacity
  for capacity in "${capacities[@]}"; do
    awk -F, -v capacity="$capacity" '
      NR > 1 && $1 == capacity {
        total++
        if ($6 == "completed") {
          completed++
          wall_sum += $11
          runtime_sum += $13
          rss_sum += $14
          budget_fallbacks += $16
          capacity_fallbacks += $17
          if ($19 > maximum_constraints)
            maximum_constraints = $19
          if ($12 ~ /^[0-9]+([.][0-9]+)?$/ && $12 > 0) {
            log_ratio_sum += log($12)
            paired++
          }
          if ($10 == "no")
            mismatch = 1
        }
      }
      END {
        wall_mean = runtime_mean = rss_mean = "NA"
        ratio = speedup = "NA"
        if (completed > 0) {
          wall_mean = sprintf("%.3f", wall_sum / completed)
          runtime_mean = sprintf("%.3f", runtime_sum / completed)
          rss_mean = sprintf("%.3f", rss_sum / completed)
        }
        if (paired > 0) {
          ratio_value = exp(log_ratio_sum / paired)
          ratio = sprintf("%.6f", ratio_value)
          speedup = sprintf("%.6f", 1 / ratio_value)
        }
        all_match = "NA"
        if (mismatch)
          all_match = "no"
        else if (completed > 0 && completed == total)
          all_match = "yes"
        printf "%s,%d,%d,%d,%s,%s,%s,%s,%s,%.0f,%.0f,%.0f,%s\n", \
          capacity, completed, total, paired, ratio, speedup, wall_mean, \
          runtime_mean, rss_mean, budget_fallbacks, capacity_fallbacks, \
          maximum_constraints, all_match
      }
    ' "$runs_csv" >>"$summary_csv"
  done
}

campaign_failed=0
for group in "${FINAL_SELECTED_GROUPS[@]}"; do
  group_root=$(final_group_result_dir "$result_root" "$group")
  mkdir -p "$group_root/logs"
  runs_csv="$group_root/results.csv"
  if [[ ! -f $runs_csv ]]; then
    printf '%s\n' \
      'capacity,graph,repetition,n,m,status,exit_code,cliques,stored_cliques,count_match_128,wall_ms,wall_ratio_vs_128,runtime_ms,max_rss_kb,ccr_states,budget_fallbacks,capacity_fallbacks,seed_solver_calls,maximum_seed_constraints,configured_capacity,dynamic,budget' \
      >"$runs_csv"
  fi

  for spec in "${FINAL_SELECTED_DATASETS[@]}"; do
    IFS='|' read -r spec_group graph n m pure_input _ <<<"$spec"
    [[ $spec_group == "$group" ]] || continue
    for ((repetition = 1; repetition <= repetitions; repetition++)); do
      for capacity in "${run_capacities[@]}"; do
        run_one "$capacity" "$graph" "$repetition" "$n" "$m" \
          "$pure_input" "$group_root" "$runs_csv"
      done
    done
  done

  write_summary "$runs_csv" "$group_root/summary.csv"
  if ! awk -F, 'NR > 1 && ($6 != "completed" || $10 == "no") { exit 1 }' \
      "$runs_csv"; then
    campaign_failed=1
  fi
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${CAPACITY_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
