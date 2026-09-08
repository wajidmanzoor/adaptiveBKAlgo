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
interval_us=${MEMORY_INTERVAL_US:-1000}
timeout_seconds=${MEMORY_TIMEOUT_SECONDS:-3600}
build_jobs=${MEMORY_BUILD_JOBS:-4}
dataset_filter=${MEMORY_DATASETS:-}
pure_budget=${MEMORY_PURE_BUDGET:-1000}
system_spec=${MEMORY_SYSTEMS:-comparison}

final_require_positive_integer MEMORY_INTERVAL_US "$interval_us" || exit 2
final_require_positive_integer MEMORY_TIMEOUT_SECONDS "$timeout_seconds" || exit 2
final_require_positive_integer MEMORY_BUILD_JOBS "$build_jobs" || exit 2
if [[ $pure_budget != unlimited ]]; then
  final_require_nonnegative_integer MEMORY_PURE_BUDGET "$pure_budget" || exit 2
fi
final_require_tools cmake awk sha256sum uname find sort || exit 2
if [[ ! -r /proc/self/status ]]; then
  echo "The memory ablation requires Linux /proc process status data." >&2
  exit 2
fi

case "$system_spec" in
  comparison|normal)
    requested_systems=(pure hbbmc_plus_plus)
    ;;
  all)
    requested_systems=(pure hbbmc hbbmc_et1 hbbmc_et2 hbbmc_et3 hbbmc_plus hbbmc_gr_et1 hbbmc_gr_et2 hbbmc_plus_plus)
    ;;
  *)
    IFS=',' read -r -a requested_systems <<<"$system_spec"
    ;;
esac

declare -A selected_systems
for system in "${requested_systems[@]}"; do
  case "$system" in
    pure|hbbmc|hbbmc_et1|hbbmc_et2|hbbmc_et3|hbbmc_plus|hbbmc_gr_et1|hbbmc_gr_et2|hbbmc_plus_plus)
      selected_systems[$system]=1
      ;;
    *)
      echo "Unknown MEMORY_SYSTEMS entry: $system" >&2
      exit 2
      ;;
  esac
done
if [[ ${#selected_systems[@]} -eq 0 ]]; then
  echo "MEMORY_SYSTEMS selected no programs." >&2
  exit 2
fi

systems=()
if [[ -n ${selected_systems[pure]:-} ]]; then
  systems+=(pure)
fi
for system in "${requested_systems[@]}"; do
  [[ $system == pure ]] && continue
  if [[ -n ${selected_systems[$system]:-} ]]; then
    systems+=("$system")
    unset 'selected_systems[$system]'
  fi
done

final_collect_grouped_datasets "$data_root" "$dataset_filter" 1 || exit 2
mkdir -p "$result_root/build_logs"
sampler_build="$result_root/build/rss_sampler"
pure_build="$result_root/build/pure"
hbbmc_build="$result_root/build/hbbmc"
pure_source="$script_dir/pure"
hbbmc_source="$script_dir/hbbmc"
final_build rss_sampler "$script_dir" "$sampler_build" \
  "$result_root/build_logs" "$build_jobs" || exit 2
final_build pure_memory "$pure_source" "$pure_build" \
  "$result_root/build_logs" "$build_jobs" || exit 2
final_build hbbmc_memory "$hbbmc_source" "$hbbmc_build" \
  "$result_root/build_logs" "$build_jobs" -DBUILD_TESTING=OFF || exit 2

sampler="$sampler_build/rss_sampler"
pure_binary="$pure_build/$FINAL_PURE_BINARY_NAME"
hbbmc_binary="$hbbmc_build/$FINAL_HBBMC_BINARY_NAME"
if [[ ! -x $sampler || ! -x $pure_binary || ! -x $hbbmc_binary ]]; then
  echo "A memory-ablation executable is missing." >&2
  exit 2
fi

{
  echo "campaign=Pure versus HBBMC++ memory traces"
  echo "data_root=$data_root"
  echo "adjacency_root=$FINAL_ADJACENCY_ROOT"
  echo "edge_root=$FINAL_EDGE_ROOT"
  echo "result_root=$result_root"
  echo "groups=${FINAL_SELECTED_GROUPS[*]}"
  echo "graphs=${#FINAL_SELECTED_DATASETS[@]}"
  echo "interval_us=$interval_us"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "execution=sequential"
  echo "systems=${systems[*]}"
  echo "pure_budget=$pure_budget"
  echo "pure_small_q_full_pxr_threshold=4"
  echo "pure_early_termination=all terminals enabled"
  echo "pure_pruning=all seven rules enabled"
  echo "hbbmc_plus_plus=RMCE graph reduction and ET level 3"
  echo "minimum_clique_size=3"
  echo "sampling_metric=VmRSS and VmHWM from /proc/PID/status"
  echo "sampling_scope=direct algorithm process including input parsing"
  echo "sampling_schedule=best-effort monotonic deadlines; elapsed_us records actual sample time"
  echo "clique_storage=unconditional in Pure and every HBBMC configuration"
  sha256sum "$sampler" "$pure_binary" "$hbbmc_binary"
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local runs_csv=$1
  local system=$2
  local graph=$3
  awk -F, -v selected="$system" -v graph="$graph" \
    'NR > 1 && $1 == selected && $2 == graph { found = 1 }
     END { exit !found }' "$runs_csv"
}

pure_reference_count() {
  local runs_csv=$1
  local graph=$2
  awk -F, -v graph="$graph" \
    'NR > 1 && $1 == "pure" && $2 == graph && $5 == "completed" {
       value = $7
     }
     END { print value }' "$runs_csv"
}

run_sampled() {
  local label=$1
  local stdout_file=$2
  local stderr_file=$3
  shift 3
  local start_ns run_pid heartbeat_pid now_ns

  echo "START $label interval=${interval_us}us timeout=${timeout_seconds}s"
  start_ns=$(date +%s%N)
  "$@" >"$stdout_file" 2>"$stderr_file" &
  run_pid=$!
  (
    while true; do
      sleep 30
      kill -0 "$run_pid" 2>/dev/null || exit
      now_ns=$(date +%s%N)
      echo "HEARTBEAT $label elapsed=$(((now_ns - start_ns) / 1000000000))s"
    done
  ) &
  heartbeat_pid=$!
  wait "$run_pid"
  MEMORY_RUN_EXIT_CODE=$?
  kill "$heartbeat_pid" 2>/dev/null || true
  wait "$heartbeat_pid" 2>/dev/null || true
}

run_one() {
  local system=$1
  local graph=$2
  local n=$3
  local m=$4
  local pure_input=$5
  local hbbmc_input=$6
  local group_root=$7
  local runs_csv=$8
  local safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  local trace_dir="$group_root/traces/$safe_graph"
  local trace_file="$trace_dir/${system}.txt"
  if already_recorded "$runs_csv" "$system" "$graph" &&
      [[ -f $trace_file ]]; then
    echo "SKIP system=$system graph=$graph reason=recorded"
    return
  fi

  local graph_reduction=none
  local et=enabled
  local expected_plus_plus=false
  local input=$pure_input
  local -a program=(env OMP_NUM_THREADS=1 "$pure_binary" "$input" --method optimized --budget "$pure_budget" --min-clique-size 3 --small-q-pxr 4)
  case "$system" in
    pure)
      ;;
    hbbmc)
      graph_reduction=none
      et=0
      ;;
    hbbmc_et1)
      graph_reduction=none
      et=1
      ;;
    hbbmc_et2)
      graph_reduction=none
      et=2
      ;;
    hbbmc_et3)
      graph_reduction=none
      et=3
      ;;
    hbbmc_plus)
      graph_reduction=rmce
      et=0
      ;;
    hbbmc_gr_et1)
      graph_reduction=rmce
      et=1
      ;;
    hbbmc_gr_et2)
      graph_reduction=rmce
      et=2
      ;;
    hbbmc_plus_plus)
      graph_reduction=rmce
      et=3
      expected_plus_plus=true
      ;;
  esac
  if [[ $system != pure ]]; then
    input=$hbbmc_input
    program=(env OMP_NUM_THREADS=1 "$hbbmc_binary" "$input" --graph-reduction "$graph_reduction" --et "$et" --num-vertices "$n" --min-clique-size 3)
  fi

  local log_dir="$group_root/logs/$safe_graph"
  local base="$log_dir/$system"
  local stdout_file="$base.stdout"
  local stderr_file="$base.stderr"
  local -a sampled_command=("$sampler" --interval-us "$interval_us" --timeout-seconds "$timeout_seconds" --output "$trace_file" -- "${program[@]}")
  mkdir -p "$log_dir" "$trace_dir"
  final_write_command "$base.command" "${sampled_command[@]}"
  run_sampled "system=$system graph=$graph" "$stdout_file" "$stderr_file" \
    "${sampled_command[@]}"

  local exit_code=$MEMORY_RUN_EXIT_CODE
  local samples sampled_peak hwm_peak wait4_peak elapsed sampler_interval
  local timed_out child_exit child_signal trace_rows
  samples=$(final_output_value "$stdout_file" sampler.samples)
  sampled_peak=$(final_output_value "$stdout_file" sampler.peak_sampled_rss_kb)
  hwm_peak=$(final_output_value "$stdout_file" sampler.peak_observed_hwm_kb)
  wait4_peak=$(final_output_value "$stdout_file" sampler.wait4_peak_rss_kb)
  elapsed=$(final_output_value "$stdout_file" sampler.elapsed_us)
  sampler_interval=$(final_output_value "$stdout_file" sampler.interval_us)
  timed_out=$(final_output_value "$stdout_file" sampler.timed_out)
  child_exit=$(final_output_value "$stdout_file" sampler.child_exit_code)
  child_signal=$(final_output_value "$stdout_file" sampler.child_signal)
  trace_rows=0
  if [[ -f $trace_file ]]; then
    trace_rows=$(awk '$1 !~ /^#/ && NF == 3 { rows++ }
                      END { print rows + 0 }' "$trace_file")
  fi

  local cliques stored runtime status algorithm_valid=no
  if [[ $system == pure ]]; then
    cliques=$(final_output_value "$stdout_file" pure.cliques)
    stored=$(final_output_value "$stdout_file" pure.stored_cliques)
    runtime=$(final_output_value "$stdout_file" pure.runtime_ms)
    if [[ $(final_output_value "$stdout_file" pure.config.budget) == "$pure_budget" &&
          $(final_output_value "$stdout_file" pure.config.small_q_full_pxr_threshold) == 4 &&
          $(final_output_value "$stdout_file" pure.config.et) == 1 &&
          $(final_output_value "$stdout_file" pure.config.et.findone_2plex) == 1 &&
          $(final_output_value "$stdout_file" pure.config.et.findone_3plex) == 1 &&
          $(final_output_value "$stdout_file" pure.config.et.full_pxr_2plex) == 1 &&
          $(final_output_value "$stdout_file" pure.config.et.full_pxr_3plex) == 1 ]] &&
        final_pruning_config_matches "$stdout_file" none; then
      algorithm_valid=yes
    fi
  else
    cliques=$(final_output_value "$stdout_file" maximal_cliques)
    stored=$(final_output_value "$stdout_file" stored_cliques)
    runtime=$(final_output_value "$stdout_file" algorithm_runtime_ms)
    if [[ $(final_output_value "$stdout_file" clique_storage) == all &&
          $(final_output_value "$stdout_file" graph_reduction) == "$graph_reduction" &&
          $(final_output_value "$stdout_file" early_termination_threshold) == "$et" &&
          $(final_output_value "$stdout_file" hbbmc_plus_plus) == "$expected_plus_plus" ]]; then
      algorithm_valid=yes
    fi
  fi

  if [[ $exit_code -eq 0 ]] &&
      final_all_uint "$cliques" "$stored" "$samples" "$sampled_peak" \
        "$hwm_peak" "$wait4_peak" "$elapsed" "$sampler_interval" \
        "$timed_out" "$child_exit" "$child_signal" "$trace_rows" &&
      final_all_number "$runtime" &&
      [[ $samples -gt 0 && $samples -eq $trace_rows &&
         $stored == "$cliques" && $sampler_interval == "$interval_us" &&
         $timed_out == 0 && $child_exit == 0 && $child_signal == 0 &&
         $algorithm_valid == yes ]]; then
    status=completed
  else
    status=$(final_status_from_exit "$exit_code")
  fi

  local reference_count count_match=NA
  if [[ $system == pure && $status == completed ]]; then
    reference_count=$cliques
    count_match=yes
  else
    reference_count=$(pure_reference_count "$runs_csv" "$graph")
    if [[ $status == completed && -n $reference_count ]]; then
      [[ $cliques == "$reference_count" ]] && count_match=yes || count_match=no
    fi
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$system" "$graph" "$n" "$m" "$status" "$exit_code" "$cliques" \
    "$stored" "$count_match" "$interval_us" "$samples" "$sampled_peak" \
    "$hwm_peak" "$wait4_peak" "$elapsed" "$runtime" \
    "$graph_reduction" "$et" "$trace_file" >>"$runs_csv"
  echo "DONE system=$system graph=$graph status=$status cliques=${cliques:-NA} stored=${stored:-NA} samples=${samples:-NA} peak_rss_kb=${sampled_peak:-NA} match=$count_match"
}

campaign_failed=0
for group in "${FINAL_SELECTED_GROUPS[@]}"; do
  group_root=$(final_group_result_dir "$result_root" "$group")
  mkdir -p "$group_root/logs" "$group_root/traces"
  runs_csv="$group_root/results.csv"
  if [[ ! -f $runs_csv ]]; then
    printf '%s\n' \
      'system,graph,n,m,status,exit_code,cliques,stored_cliques,count_match_pure,interval_us,samples,peak_sampled_rss_kb,peak_observed_hwm_kb,wait4_peak_rss_kb,elapsed_us,algorithm_runtime_ms,graph_reduction,et,trace_file' \
      >"$runs_csv"
  fi

  for system in "${systems[@]}"; do
    for spec in "${FINAL_SELECTED_DATASETS[@]}"; do
      IFS='|' read -r spec_group graph n m pure_input hbbmc_input <<<"$spec"
      [[ $spec_group == "$group" ]] || continue
      run_one "$system" "$graph" "$n" "$m" "$pure_input" \
        "$hbbmc_input" "$group_root" "$runs_csv"
    done
  done

  if ! awk -F, 'NR > 1 && ($5 != "completed" || $9 == "no") { exit 1 }' \
      "$runs_csv"; then
    campaign_failed=1
  fi
done

echo "CAMPAIGN_COMPLETE results=$result_root"
if [[ ${MEMORY_FAIL_ON_ERROR:-0} == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
