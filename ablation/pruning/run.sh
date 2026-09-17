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
timeout_seconds=${PRUNING_TIMEOUT_SECONDS:-3600}
build_jobs=${PRUNING_BUILD_JOBS:-4}
dataset_filter=${PRUNING_DATASETS:-}
budget=${PRUNING_BUDGET:-1000}

final_require_positive_integer PRUNING_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer PRUNING_BUILD_JOBS "$build_jobs" || exit 2
if [[ $budget != unlimited ]]; then
  final_require_nonnegative_integer PRUNING_BUDGET "$budget" || exit 2
fi
final_require_tools cmake timeout /usr/bin/time awk sha256sum uname find sort || exit 2
final_collect_grouped_datasets "$data_root" "$dataset_filter" 0 || exit 2

mkdir -p "$result_root/build_logs"
declare -A binaries
variants=(all_rules)
pure_source="$script_dir/pure"
for rule in "${FINAL_PRUNING_RULES[@]}"; do
  variants+=("leave_$rule")
done

for variant in "${variants[@]}"; do
  build_root="$result_root/build/$variant"
  if [[ $variant == all_rules ]]; then
    final_build "$variant" "$pure_source" "$build_root" \
      "$result_root/build_logs" "$build_jobs" || exit 2
  else
    rule=${variant#leave_}
    option="-DPURE_DISABLE_PRUNING_${rule^^}=ON"
    final_build "$variant" "$pure_source" "$build_root" \
      "$result_root/build_logs" "$build_jobs" "$option" || exit 2
  fi
  binaries[$variant]="$build_root/$FINAL_PURE_BINARY_NAME"
  [[ -x ${binaries[$variant]} ]] || {
    echo "Missing reorder binary for $variant." >&2
    exit 2
  }
done

{
  echo "campaign=final Reorder+CCRMCE seven-rule leave-one-out pruning ablation"
  echo "data_root=$data_root"
  echo "adjacency_root=$FINAL_ADJACENCY_ROOT"
  echo "result_root=$result_root"
  echo "groups=${FINAL_SELECTED_GROUPS[*]}"
  echo "graphs=${#FINAL_SELECTED_DATASETS[@]}"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "execution=sequential"
  echo "minimum_clique_size=3"
  echo "small_q_ccr_threshold=32"
  echo "budget=$budget"
  echo "hitset_capacity=128"
  echo "early_termination=ET1 ET2 ET3 disabled"
  echo "rules=${FINAL_PRUNING_RULES[*]}"
  echo "baseline=all rules enabled"
  for variant in "${variants[@]}"; do
    sha256sum "${binaries[$variant]}"
  done
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local runs_csv=$1
  local variant=$2
  local graph=$3
  awk -F, -v variant="$variant" -v graph="$graph" \
    'NR > 1 && $1 == variant && $3 == graph { found = 1 }
     END { exit !found }' "$runs_csv"
}

baseline_field() {
  local runs_csv=$1
  local graph=$2
  local field=$3
  awk -F, -v graph="$graph" -v field="$field" \
    'NR > 1 && $1 == "all_rules" && $3 == graph { value = $field }
     END { print value }' "$runs_csv"
}

run_one() {
  local variant=$1
  local disabled=$2
  local graph=$3
  local n=$4
  local m=$5
  local input=$6
  local group_root=$7
  local runs_csv=$8
  if already_recorded "$runs_csv" "$variant" "$graph"; then
    echo "SKIP variant=$variant graph=$graph reason=recorded"
    return
  fi

  local binary=${binaries[$variant]}
  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local log_dir="$group_root/logs/$safe_graph"
  local base="$log_dir/$variant"
  mkdir -p "$log_dir"
  local -a command=(
    env OMP_NUM_THREADS=1 "$binary" "$input"
    --budget "$budget" --min-clique-size 3
  )
  final_write_command "$base.command" "${command[@]}"
  final_run_timed "variant=$variant graph=$graph" "$base.stdout" \
    "$base.stderr" "$base.resources" "$timeout_seconds" "${command[@]}"

  local exit_code=$FINAL_RUN_EXIT_CODE
  local wall_ms=$FINAL_RUN_WALL_MS
  local cliques runtime_ms findone_states full_states ccr_states fallbacks configured status
  local reference_count reference_wall count_match=NA slowdown=NA
  cliques=$(final_output_value "$base.stdout" reorder.cliques)
  runtime_ms=$(final_output_value "$base.stdout" reorder.runtime_ms)
  findone_states=$(final_output_value "$base.stdout" reorder.ccr.findone_states)
  full_states=$(final_output_value "$base.stdout" reorder.ccr.full_states)
  ccr_states=
  if final_all_uint "$findone_states" "$full_states"; then
    ccr_states=$((findone_states + full_states))
  fi
  fallbacks=$(final_output_value "$base.stdout" reorder.budget_fallbacks)
  configured=$(final_output_value "$base.stdout" reorder.budget)

  if [[ $exit_code -eq 0 ]] &&
      final_all_uint "$cliques" "$ccr_states" "$fallbacks" &&
      final_all_number "$runtime_ms" &&
      [[ $configured == "$budget" ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.small_q_ccr_threshold) == 32 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et1) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et2) == 0 ]] &&
      [[ $(final_output_value "$base.stdout" reorder.config.et3) == 0 ]] &&
      final_pruning_config_matches "$base.stdout" "$disabled"; then
    status=completed
  else
    status=$(final_status_from_exit "$exit_code")
  fi

  if [[ $variant == all_rules ]]; then
    reference_count=$cliques
    reference_wall=$wall_ms
  else
    reference_count=$(baseline_field "$runs_csv" "$graph" 8)
    reference_wall=$(baseline_field "$runs_csv" "$graph" 9)
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
    "$variant" "$disabled" "$graph" "$n" "$m" "$status" \
    "$exit_code" "$cliques" "$wall_ms" "$runtime_ms" "$ccr_states" \
    "$fallbacks" "$reference_count" "$count_match" "$slowdown" \
    >>"$runs_csv"
  echo "DONE variant=$variant graph=$graph status=$status cliques=${cliques:-NA} match=$count_match"
}

campaign_failed=0
for group in "${FINAL_SELECTED_GROUPS[@]}"; do
  group_root=$(final_group_result_dir "$result_root" "$group")
  mkdir -p "$group_root/logs"
  runs_csv="$group_root/results.csv"
  if [[ ! -f $runs_csv ]]; then
    printf '%s\n' \
      'variant,disabled_rule,graph,n,m,status,exit_code,cliques,wall_ms,runtime_ms,ccr_states,budget_fallbacks,reference_count,count_match,slowdown_vs_all_rules' \
      >"$runs_csv"
  fi

  for variant in "${variants[@]}"; do
    disabled=none
    [[ $variant != all_rules ]] && disabled=${variant#leave_}
    for spec in "${FINAL_SELECTED_DATASETS[@]}"; do
      IFS='|' read -r spec_group graph n m pure_input _ <<<"$spec"
      [[ $spec_group == "$group" ]] || continue
      run_one "$variant" "$disabled" "$graph" "$n" "$m" "$pure_input" \
        "$group_root" "$runs_csv"
    done
  done

  if ! awk -F, 'NR > 1 && ($6 != "completed" || $14 == "no") { exit 1 }' \
      "$runs_csv"; then
    campaign_failed=1
  fi
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${PRUNING_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
