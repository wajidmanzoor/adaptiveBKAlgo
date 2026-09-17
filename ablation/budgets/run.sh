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
timeout_seconds=${BUDGET_TIMEOUT_SECONDS:-3600}
build_jobs=${BUDGET_BUILD_JOBS:-4}
dataset_filter=${BUDGET_DATASETS:-}
reference_budget=${BUDGET_REFERENCE:-10000}
budget_csv=${BUDGET_VALUES:-0,100,1000,10000,100000,unlimited}

final_require_positive_integer BUDGET_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer BUDGET_BUILD_JOBS "$build_jobs" || exit 2
final_require_tools cmake timeout /usr/bin/time awk sha256sum uname find sort || exit 2

IFS=',' read -r -a requested_budgets <<<"$budget_csv"
budgets=("$reference_budget")
reference_seen=0
for budget in "${requested_budgets[@]}"; do
  [[ -n $budget ]] || {
    echo "BUDGET_VALUES contains an empty value." >&2
    exit 2
  }
  if [[ $budget != unlimited ]]; then
    final_require_nonnegative_integer budget "$budget" || exit 2
  fi
  [[ $budget == "$reference_budget" ]] && reference_seen=1
  [[ $budget == "$reference_budget" ]] && continue
  budgets+=("$budget")
done
if [[ $reference_seen -ne 1 ]]; then
  echo "BUDGET_REFERENCE must also appear in BUDGET_VALUES." >&2
  exit 2
fi

final_collect_grouped_datasets "$data_root" "$dataset_filter" 0 || exit 2
mkdir -p "$result_root/build_logs"
build_root="$result_root/build/pure"
pure_source="$script_dir/pure"
final_build pure_budget "$pure_source" "$build_root" \
  "$result_root/build_logs" "$build_jobs" || exit 2
binary="$build_root/$FINAL_PURE_BINARY_NAME"
[[ -x $binary ]] || {
  echo "Missing reorder binary: $binary" >&2
  exit 2
}

{
  echo "campaign=final Reorder+CCRMCE budget sweep"
  echo "data_root=$data_root"
  echo "adjacency_root=$FINAL_ADJACENCY_ROOT"
  echo "result_root=$result_root"
  echo "groups=${FINAL_SELECTED_GROUPS[*]}"
  echo "graphs=${#FINAL_SELECTED_DATASETS[@]}"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "execution=sequential"
  echo "minimum_clique_size=3"
  echo "small_q_ccr_threshold=32"
  echo "adaptive_direct_q_threshold=256"
  echo "hitset_capacity=128"
  echo "early_termination=ET1 ET2 ET3 disabled"
  echo "pruning=production profile (subsumption disabled; six other rules enabled)"
  echo "reference_budget=$reference_budget"
  echo "budgets=${budgets[*]}"
  echo "budget_semantics=numeric N permits N compatibility examinations per seed-solver call; 0 forces fallback at the first examination; unlimited imposes no limit"
  sha256sum "$binary"
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local runs_csv=$1
  local budget=$2
  local graph=$3
  awk -F, -v budget="$budget" -v graph="$graph" \
    'NR > 1 && $1 == budget && $2 == graph { found = 1 }
     END { exit !found }' "$runs_csv"
}

csv_field() {
  local runs_csv=$1
  local budget=$2
  local graph=$3
  local field=$4
  awk -F, -v budget="$budget" -v graph="$graph" -v field="$field" \
    'NR > 1 && $1 == budget && $2 == graph { value = $field }
     END { print value }' "$runs_csv"
}

run_one() {
  local budget=$1
  local graph=$2
  local n=$3
  local m=$4
  local input=$5
  local group_root=$6
  local runs_csv=$7
  if already_recorded "$runs_csv" "$budget" "$graph"; then
    echo "SKIP budget=$budget graph=$graph reason=recorded"
    return
  fi

  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local log_dir="$group_root/logs/$safe_graph"
  local base="$log_dir/budget_$budget"
  mkdir -p "$log_dir"
  local -a command=(
    env OMP_NUM_THREADS=1 "$binary" "$input"
    --budget "$budget" --min-clique-size 3
  )
  final_write_command "$base.command" "${command[@]}"
  final_run_timed "budget=$budget graph=$graph" "$base.stdout" \
    "$base.stderr" "$base.resources" "$timeout_seconds" "${command[@]}"

  local exit_code=$FINAL_RUN_EXIT_CODE
  local wall_ms=$FINAL_RUN_WALL_MS
  local cliques runtime_ms findone_states full_states ccr_states
  local fallbacks capacity_fallbacks configured status
  local reference_count reference_wall count_match=NA slowdown=NA
  cliques=$(final_output_value "$base.stdout" reorder.cliques)
  runtime_ms=$(final_output_value "$base.stdout" reorder.runtime_ms)
  findone_states=$(final_output_value "$base.stdout" reorder.ccr.findone_states)
  full_states=$(final_output_value "$base.stdout" reorder.ccr.full_states)
  fallbacks=$(final_output_value "$base.stdout" reorder.budget_fallbacks)
  capacity_fallbacks=$(final_output_value "$base.stdout" reorder.capacity_fallbacks)
  configured=$(final_output_value "$base.stdout" reorder.budget)
  ccr_states=
  if final_all_uint "$findone_states" "$full_states"; then
    ccr_states=$((findone_states + full_states))
  fi

  if [[ $exit_code -eq 0 ]] &&
      final_all_uint "$cliques" "$ccr_states" "$fallbacks" "$capacity_fallbacks" &&
      final_all_number "$runtime_ms" &&
      [[ $configured == "$budget" && $capacity_fallbacks == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.small_q_ccr_threshold) == 32 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et1) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et2) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et3) == 0 ]] &&
      final_pruning_config_matches "$base.stdout" production; then
    status=completed
  else
    status=$(final_status_from_exit "$exit_code")
  fi

  if [[ $budget == "$reference_budget" ]]; then
    reference_count=$cliques
    reference_wall=$wall_ms
  else
    reference_count=$(csv_field "$runs_csv" "$reference_budget" "$graph" 7)
    reference_wall=$(csv_field "$runs_csv" "$reference_budget" "$graph" 8)
  fi
  if [[ $status == completed && -n $reference_count ]]; then
    [[ $cliques == "$reference_count" ]] && count_match=yes || count_match=no
  fi
  if [[ $status == completed && -n $reference_wall &&
        $reference_wall =~ ^[0-9]+$ && $reference_wall -gt 0 ]]; then
    slowdown=$(awk -v wall="$wall_ms" -v base="$reference_wall" \
      'BEGIN { printf "%.3f", wall / base }')
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$budget" "$graph" "$n" "$m" "$status" "$exit_code" "$cliques" \
    "$wall_ms" "$runtime_ms" "$ccr_states" "$fallbacks" "$configured" \
    "$reference_count" "$count_match" "$slowdown" >>"$runs_csv"
  echo "DONE budget=$budget graph=$graph status=$status cliques=${cliques:-NA} match=$count_match fallbacks=${fallbacks:-NA}"
}

campaign_failed=0
for group in "${FINAL_SELECTED_GROUPS[@]}"; do
  group_root=$(final_group_result_dir "$result_root" "$group")
  mkdir -p "$group_root/logs"
  runs_csv="$group_root/results.csv"
  if [[ ! -f $runs_csv ]]; then
    printf '%s\n' \
      'budget,graph,n,m,status,exit_code,cliques,wall_ms,runtime_ms,ccr_states,budget_fallbacks,configured_budget,reference_count,count_match,slowdown_vs_reference' \
      >"$runs_csv"
  fi

  for budget in "${budgets[@]}"; do
    for spec in "${FINAL_SELECTED_DATASETS[@]}"; do
      IFS='|' read -r spec_group graph n m pure_input _ <<<"$spec"
      [[ $spec_group == "$group" ]] || continue
      run_one "$budget" "$graph" "$n" "$m" "$pure_input" \
        "$group_root" "$runs_csv"
    done
  done

  if ! awk -F, 'NR > 1 && ($5 != "completed" || $14 == "no") { exit 1 }' \
      "$runs_csv"; then
    campaign_failed=1
  fi
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${BUDGET_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
