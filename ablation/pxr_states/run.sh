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
timeout_seconds=${PXR_TIMEOUT_SECONDS:-3600}
build_jobs=${PXR_BUILD_JOBS:-4}
dataset_filter=${PXR_DATASETS:-}
reorder_budget=${PXR_REORDER_BUDGET:-1000}

final_require_positive_integer PXR_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer PXR_BUILD_JOBS "$build_jobs" || exit 2
if [[ $reorder_budget != unlimited ]]; then
  final_require_nonnegative_integer PXR_REORDER_BUDGET "$reorder_budget" || exit 2
fi
final_require_tools cmake timeout /usr/bin/time awk sha256sum uname find sort || exit 2
final_collect_grouped_datasets "$data_root" "$dataset_filter" 1 || exit 2

mkdir -p "$result_root/build_logs"
pure_build="$result_root/build/reorder_no_et"
hbbmc_build="$result_root/build/hbbmc_no_reduction_no_et"
comparison_build="$result_root/build/tomita_ccrmce"
pure_source="$script_dir/pure"
hbbmc_source="$script_dir/hbbmc"
final_build reorder_no_et "$pure_source" "$pure_build" \
  "$result_root/build_logs" "$build_jobs" || exit 2
final_build hbbmc_no_reduction_no_et "$hbbmc_source" "$hbbmc_build" \
  "$result_root/build_logs" "$build_jobs" -DBUILD_TESTING=OFF || exit 2
final_build tomita_ccrmce_states "$script_dir" "$comparison_build" \
  "$result_root/build_logs" "$build_jobs" || exit 2
pure_binary="$pure_build/$FINAL_PURE_BINARY_NAME"
hbbmc_binary="$hbbmc_build/$FINAL_HBBMC_BINARY_NAME"
tomita_adapter="$comparison_build/tomita_stream_adapter"
tomita_binary="$comparison_build/tomita_states"
ccrmce_binary="$comparison_build/ccrmce_states"
if [[ ! -x $pure_binary || ! -x $hbbmc_binary || ! -x $tomita_adapter ||
      ! -x $tomita_binary || ! -x $ccrmce_binary ]]; then
  echo "A recursive-state executable is missing." >&2
  exit 2
fi

{
  echo "campaign=four-system exact maximal-clique recursive states"
  echo "data_root=$data_root"
  echo "adjacency_root=$FINAL_ADJACENCY_ROOT"
  echo "edge_root=$FINAL_EDGE_ROOT"
  echo "result_root=$result_root"
  echo "groups=${FINAL_SELECTED_GROUPS[*]}"
  echo "graphs=${#FINAL_SELECTED_DATASETS[@]}"
  echo "timeout_seconds_per_program=$timeout_seconds"
  echo "execution=sequential"
  echo "minimum_clique_size=3"
  echo "reorder_budget=$reorder_budget"
  echo "reorder_small_q_ccr_threshold=32"
  echo "reorder_hitset_capacity=128"
  echo "reorder_pruning=production profile"
  echo "reorder_legacy_early_termination=disabled"
  echo "hbbmc_graph_reduction=none"
  echo "hbbmc_early_termination=0"
  echo "state_definition.reorder=one recursive CCRMCE entry in FindOne or exhaustive enumeration"
  echo "state_definition.hbbmc=counter.vertex_recursive_calls"
  echo "state_definition.tomita=one entry to listAllMaximalCliquesAdjacencyListRecursive"
  echo "state_definition.ccrmce=one entry to BKTomitaRecCCRV3"
  echo "state_scope=implementation-native recursive search entries; preprocessing and root setup excluded"
  echo "clique_storage=all four systems retain every reported vertex list"
  sha256sum "$pure_binary" "$hbbmc_binary" "$tomita_adapter" \
    "$tomita_binary" "$ccrmce_binary"
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local runs_csv=$1
  local graph=$2
  awk -F, -v graph="$graph" \
    'NR > 1 && $1 == graph { found = 1 } END { exit !found }' "$runs_csv"
}
final_colon_value() {
  local file=$1
  local key=$2
  awk -F: -v key="$key" '$1 == key {
      value = substr($0, length($1) + 2)
    }
    END { print value }' "$file"
}


run_dataset() {
  local graph=$1
  local n=$2
  local m=$3
  local pure_input=$4
  local hbbmc_input=$5
  local group_root=$6
  local runs_csv=$7
  if already_recorded "$runs_csv" "$graph"; then
    echo "SKIP graph=$graph reason=recorded"
    return
  fi

  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local log_dir="$group_root/logs/$safe_graph"
  local pure_base="$log_dir/reorder"
  local hbbmc_base="$log_dir/hbbmc"
  local tomita_base="$log_dir/tomita"
  local ccrmce_base="$log_dir/ccrmce"
  local tomita_input="$log_dir/tomita.input"
  mkdir -p "$log_dir"

  final_write_command "$log_dir/tomita_adapter.command" \
    "$tomita_adapter" "$hbbmc_input" "$n" "$m"
  if ! "$tomita_adapter" "$hbbmc_input" "$n" "$m" \
      >"$tomita_input" 2>"$log_dir/tomita_adapter.stderr"; then
    echo "Tomita input conversion failed for $graph." >&2
    exit 2
  fi

  local -a pure_command=(
    env OMP_NUM_THREADS=1 "$pure_binary" "$pure_input"
    --budget "$reorder_budget" --min-clique-size 3
  )
  local -a hbbmc_command=(
    env OMP_NUM_THREADS=1 "$hbbmc_binary" "$hbbmc_input"
    --graph-reduction none --et 0 --num-vertices "$n"
    --min-clique-size 3 --counters
  )
  local -a tomita_command=(
    env OMP_NUM_THREADS=1 "$tomita_binary" "$tomita_input"
  )
  local -a ccrmce_command=(
    env OMP_NUM_THREADS=1 "$ccrmce_binary" noUVM -f_txt "$hbbmc_input"
  )

  final_write_command "$pure_base.command" "${pure_command[@]}"
  final_run_timed "system=reorder graph=$graph" "$pure_base.stdout" \
    "$pure_base.stderr" "$pure_base.resources" "$timeout_seconds" \
    "${pure_command[@]}"
  local pure_exit=$FINAL_RUN_EXIT_CODE
  local pure_wall=$FINAL_RUN_WALL_MS

  local reorder_cliques reorder_stored pure_pxr pure_findone pure_full pure_et pure_total
  local pure_runtime pure_checks pure_fallbacks reorder_status reorder_invariant=no
  reorder_cliques=$(final_output_value "$pure_base.stdout" reorder.cliques)
  reorder_stored=$(final_output_value "$pure_base.stdout" reorder.stored_cliques)
  pure_pxr=
  pure_findone=$(final_output_value "$pure_base.stdout" reorder.ccr.findone_states)
  pure_full=$(final_output_value "$pure_base.stdout" reorder.ccr.full_states)
  pure_et=0
  pure_total=
  if final_all_uint "$pure_findone" "$pure_full"; then
    pure_pxr=$((pure_findone + pure_full))
    pure_total=$pure_pxr
  fi
  pure_runtime=$(final_output_value "$pure_base.stdout" reorder.runtime_ms)
  pure_checks=$pure_pxr
  pure_fallbacks=$(final_output_value "$pure_base.stdout" reorder.budget_fallbacks)

  if [[ $pure_exit -eq 0 ]] &&
      final_all_uint "$reorder_cliques" "$reorder_stored" "$pure_pxr" "$pure_findone" \
        "$pure_full" "$pure_et" "$pure_total" "$pure_checks" \
        "$pure_fallbacks" &&
      final_all_number "$pure_runtime" &&
      [[ $reorder_stored == "$reorder_cliques" ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.minimum_clique_size) == 3 ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.budget) == "$reorder_budget" ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.config.small_q_ccr_threshold) == 32 ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.config.et1) == 0 ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.config.et2) == 0 ]] &&
      [[ $(final_output_value "$pure_base.stdout" reorder.config.et3) == 0 ]] &&
      final_pruning_config_matches "$pure_base.stdout" production; then
    if [[ $pure_pxr -eq $((pure_findone + pure_full)) &&
          $pure_checks -eq $pure_pxr && $pure_et -eq 0 &&
          $pure_total -eq $pure_pxr ]]; then
      reorder_invariant=yes
      reorder_status=completed
    else
      reorder_status=failed
    fi
  else
    reorder_status=$(final_status_from_exit "$pure_exit")
  fi
  echo "DONE system=reorder graph=$graph status=$reorder_status ccr_states=${pure_pxr:-NA}"

  final_write_command "$hbbmc_base.command" "${hbbmc_command[@]}"
  final_run_timed "system=hbbmc graph=$graph" "$hbbmc_base.stdout" \
    "$hbbmc_base.stderr" "$hbbmc_base.resources" "$timeout_seconds" \
    "${hbbmc_command[@]}"
  local hbbmc_exit=$FINAL_RUN_EXIT_CODE
  local hbbmc_wall=$FINAL_RUN_WALL_MS

  local hbbmc_cliques hbbmc_stored hbbmc_pxr hbbmc_et1 hbbmc_et2 hbbmc_et3
  local hbbmc_et hbbmc_runtime hbbmc_status hbbmc_invariant=no
  hbbmc_cliques=$(final_output_value "$hbbmc_base.stdout" maximal_cliques)
  hbbmc_stored=$(final_output_value "$hbbmc_base.stdout" stored_cliques)
  hbbmc_pxr=$(final_output_value "$hbbmc_base.stdout" counter.vertex_recursive_calls)
  hbbmc_et1=$(final_output_value "$hbbmc_base.stdout" counter.et1_outputs)
  hbbmc_et2=$(final_output_value "$hbbmc_base.stdout" counter.et2_outputs)
  hbbmc_et3=$(final_output_value "$hbbmc_base.stdout" counter.et3_outputs)
  hbbmc_runtime=$(final_output_value "$hbbmc_base.stdout" algorithm_runtime_ms)
  hbbmc_et=0
  if final_all_uint "$hbbmc_et1" "$hbbmc_et2" "$hbbmc_et3"; then
    hbbmc_et=$((hbbmc_et1 + hbbmc_et2 + hbbmc_et3))
  fi

  if [[ $hbbmc_exit -eq 0 ]] &&
      final_all_uint "$hbbmc_cliques" "$hbbmc_stored" "$hbbmc_pxr" \
        "$hbbmc_et1" "$hbbmc_et2" "$hbbmc_et3" &&
      final_all_number "$hbbmc_runtime" &&
      [[ $hbbmc_stored == "$hbbmc_cliques" ]] &&
      [[ $(final_output_value "$hbbmc_base.stdout" clique_storage) == all ]] &&
      [[ $(final_output_value "$hbbmc_base.stdout" graph_reduction) == none ]] &&
      [[ $(final_output_value "$hbbmc_base.stdout" early_termination_threshold) == 0 ]] &&
      [[ $(final_output_value "$hbbmc_base.stdout" hbbmc_plus_plus) == false ]] &&
      [[ $(final_output_value "$hbbmc_base.stdout" counter.et_checks) == 0 ]] &&
      [[ $hbbmc_et -eq 0 ]]; then
    hbbmc_invariant=yes
    hbbmc_status=completed
  else
    hbbmc_status=$(final_status_from_exit "$hbbmc_exit")
  fi
  echo "DONE system=hbbmc graph=$graph status=$hbbmc_status pxr=${hbbmc_pxr:-NA}"

  final_write_command "$tomita_base.command" "${tomita_command[@]}"
  final_run_timed "system=tomita graph=$graph" "$tomita_base.stdout" \
    "$tomita_base.stderr" "$tomita_base.resources" "$timeout_seconds" \
    "${tomita_command[@]}"
  local tomita_exit=$FINAL_RUN_EXIT_CODE
  local tomita_wall=$FINAL_RUN_WALL_MS

  local tomita_cliques tomita_stored tomita_states tomita_runtime
  local tomita_status tomita_invariant=no
  tomita_cliques=$(final_output_value "$tomita_base.stdout" maximal_cliques)
  tomita_stored=$(final_output_value "$tomita_base.stdout" stored_cliques)
  tomita_states=$(final_output_value "$tomita_base.stdout" recursive_states)
  tomita_runtime=$(final_output_value "$tomita_base.stdout" algorithm_wall_ms)
  if [[ $tomita_exit -eq 0 ]] &&
      final_all_uint "$tomita_cliques" "$tomita_stored" "$tomita_states" &&
      final_all_number "$tomita_runtime" &&
      [[ $tomita_stored == "$tomita_cliques" ]] &&
      [[ $(final_output_value "$tomita_base.stdout" algorithm) == tomita-adjacency-list ]] &&
      [[ $(final_output_value "$tomita_base.stdout" minimum_clique_size) == 3 ]] &&
      [[ $(final_output_value "$tomita_base.stdout" clique_storage) == all ]]; then
    tomita_invariant=yes
    tomita_status=completed
  else
    tomita_status=$(final_status_from_exit "$tomita_exit")
  fi
  echo "DONE system=tomita graph=$graph status=$tomita_status states=${tomita_states:-NA}"

  final_write_command "$ccrmce_base.command" "${ccrmce_command[@]}"
  final_run_timed "system=ccrmce graph=$graph" "$ccrmce_base.stdout" \
    "$ccrmce_base.stderr" "$ccrmce_base.resources" "$timeout_seconds" \
    "${ccrmce_command[@]}"
  local ccrmce_exit=$FINAL_RUN_EXIT_CODE
  local ccrmce_wall=$FINAL_RUN_WALL_MS

  local ccrmce_cliques ccrmce_stored ccrmce_states ccrmce_runtime
  local ccrmce_status ccrmce_invariant=no
  ccrmce_cliques=$(final_colon_value "$ccrmce_base.stdout" Mclique)
  ccrmce_stored=$(final_colon_value "$ccrmce_base.stdout" stored_cliques)
  ccrmce_states=$(final_colon_value "$ccrmce_base.stdout" numOfDFSSearchNodes)
  ccrmce_runtime=$(final_colon_value "$ccrmce_base.stdout" time)
  ccrmce_runtime=${ccrmce_runtime%ms}
  if [[ $ccrmce_exit -eq 0 ]] &&
      final_all_uint "$ccrmce_cliques" "$ccrmce_stored" "$ccrmce_states" &&
      final_all_number "$ccrmce_runtime" &&
      [[ $ccrmce_stored == "$ccrmce_cliques" ]] &&
      [[ $(final_colon_value "$ccrmce_base.stdout" clique_storage) == all ]]; then
    ccrmce_invariant=yes
    ccrmce_status=completed
  else
    ccrmce_status=$(final_status_from_exit "$ccrmce_exit")
  fi
  echo "DONE system=ccrmce graph=$graph status=$ccrmce_status states=${ccrmce_states:-NA}"

  local hbbmc_count_match=NA tomita_count_match=NA ccrmce_count_match=NA
  local hbbmc_ratio=NA tomita_ratio=NA ccrmce_ratio=NA all_counts_match=NA
  if [[ $reorder_status == completed && $hbbmc_status == completed ]]; then
    [[ $reorder_cliques == "$hbbmc_cliques" ]] && hbbmc_count_match=yes || hbbmc_count_match=no
    if [[ $pure_pxr -gt 0 ]]; then
      hbbmc_ratio=$(awk -v states="$hbbmc_pxr" -v base="$pure_pxr" \
        'BEGIN { printf "%.6f", states / base }')
    fi
  fi
  if [[ $reorder_status == completed && $tomita_status == completed ]]; then
    [[ $reorder_cliques == "$tomita_cliques" ]] && tomita_count_match=yes || tomita_count_match=no
    if [[ $pure_pxr -gt 0 ]]; then
      tomita_ratio=$(awk -v states="$tomita_states" -v base="$pure_pxr" \
        'BEGIN { printf "%.6f", states / base }')
    fi
  fi
  if [[ $reorder_status == completed && $ccrmce_status == completed ]]; then
    [[ $reorder_cliques == "$ccrmce_cliques" ]] && ccrmce_count_match=yes || ccrmce_count_match=no
    if [[ $pure_pxr -gt 0 ]]; then
      ccrmce_ratio=$(awk -v states="$ccrmce_states" -v base="$pure_pxr" \
        'BEGIN { printf "%.6f", states / base }')
    fi
  fi
  if [[ $reorder_status == completed && $hbbmc_status == completed &&
        $tomita_status == completed && $ccrmce_status == completed ]]; then
    all_counts_match=no
    if [[ $hbbmc_count_match == yes && $tomita_count_match == yes &&
          $ccrmce_count_match == yes ]]; then
      all_counts_match=yes
    fi
  fi

  local -a row=(
    "$graph" "$n" "$m"
    "$reorder_status" "$reorder_cliques" "$reorder_stored" "$pure_pxr"
    "$pure_findone" "$pure_full" "$pure_et" "$pure_runtime" "$pure_wall"
    "$pure_fallbacks" "$reorder_invariant"
    "$hbbmc_status" "$hbbmc_cliques" "$hbbmc_stored" "$hbbmc_pxr"
    "$hbbmc_et" "$hbbmc_runtime" "$hbbmc_wall" "$hbbmc_invariant"
    "$hbbmc_count_match" "$hbbmc_ratio"
    "$tomita_status" "$tomita_cliques" "$tomita_stored" "$tomita_states"
    "$tomita_runtime" "$tomita_wall" "$tomita_invariant"
    "$tomita_count_match" "$tomita_ratio"
    "$ccrmce_status" "$ccrmce_cliques" "$ccrmce_stored" "$ccrmce_states"
    "$ccrmce_runtime" "$ccrmce_wall" "$ccrmce_invariant"
    "$ccrmce_count_match" "$ccrmce_ratio" "$all_counts_match"
  )
  {
    printf '%s' "${row[0]}"
    printf ',%s' "${row[@]:1}"
    printf '\n'
  } >>"$runs_csv"
  echo "COMPARE graph=$graph all_counts_match=$all_counts_match hbbmc_ratio=$hbbmc_ratio tomita_ratio=$tomita_ratio ccrmce_ratio=$ccrmce_ratio"
}

campaign_failed=0
for group in "${FINAL_SELECTED_GROUPS[@]}"; do
  group_root=$(final_group_result_dir "$result_root" "$group")
  mkdir -p "$group_root/logs"
  runs_csv="$group_root/results.csv"
  if [[ ! -f $runs_csv ]]; then
    printf '%s\n' \
      'graph,n,m,reorder_status,reorder_cliques,reorder_stored_cliques,reorder_ccr_states,reorder_findone_ccr_states,reorder_full_ccr_states,reorder_legacy_et_states,reorder_runtime_ms,reorder_wall_ms,reorder_budget_fallbacks,reorder_invariant,hbbmc_status,hbbmc_cliques,hbbmc_stored_cliques,hbbmc_pxr_states,hbbmc_et_states,hbbmc_runtime_ms,hbbmc_wall_ms,hbbmc_invariant,hbbmc_count_match,hbbmc_over_reorder_ccr_ratio,tomita_status,tomita_cliques,tomita_stored_cliques,tomita_recursive_states,tomita_runtime_ms,tomita_wall_ms,tomita_invariant,tomita_count_match,tomita_over_reorder_ccr_ratio,ccrmce_status,ccrmce_cliques,ccrmce_stored_cliques,ccrmce_recursive_states,ccrmce_runtime_ms,ccrmce_wall_ms,ccrmce_invariant,ccrmce_count_match,ccrmce_over_reorder_ccr_ratio,all_counts_match' \
      >"$runs_csv"
  fi

  for spec in "${FINAL_SELECTED_DATASETS[@]}"; do
    IFS='|' read -r spec_group graph n m pure_input hbbmc_input <<<"$spec"
    [[ $spec_group == "$group" ]] || continue
    run_dataset "$graph" "$n" "$m" "$pure_input" "$hbbmc_input" \
      "$group_root" "$runs_csv"
  done

  if ! awk -F, 'NR > 1 && ($4 != "completed" || $15 != "completed" ||
                            $25 != "completed" || $34 != "completed" ||
                            $43 != "yes") { exit 1 }' "$runs_csv"; then
    campaign_failed=1
  fi
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${PXR_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
