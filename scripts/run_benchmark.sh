#!/usr/bin/env bash

set -uo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/run_benchmark.sh ADJACENCY_DIR EDGE_DIR [RESULT_ROOT]

Build and compare complete process wall time for:
  - Reorder/AdaptiveBK (budget 1000);
  - the independent paper-faithful HBBMC++ implementation (RMCE + ET3);
  - sparse adjacency-list Tomita; and
  - standalone CCRMCE.

ADJACENCY_DIR and EDGE_DIR must contain the same relative GROUP/GRAPH paths.
Adjacency inputs use the Reorder format (header "n m", then n adjacency rows).
Edge inputs are normalized headerless "u v" rows with 0 <= u < v < n.

The primary measurement is wall_time_us around the entire program invocation:
process startup, native-format parsing, preprocessing, enumeration, clique
storage performed by that implementation, output, and process teardown.
Tomita's one-time native-format conversion is prepared and logged separately.

Environment:
  COMPARE_TIMEOUT_SECONDS  Per-program timeout (default: 1800)
  COMPARE_BUILD_JOBS       Parallel build jobs (default: 4)
  COMPARE_DATASETS         Comma/space-separated group, graph, stem, or path
  COMPARE_FAIL_ON_ERROR    Return nonzero after the campaign if any run failed
                           or counts disagreed (default: 0)

Each failed or timed-out graph/program row is recorded; later programs and
graphs always continue. Reuse RESULT_ROOT to resume recorded rows.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  usage
  exit 0
fi
if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
adjacency_dir=$1
edge_dir=$2
timestamp=$(date +%Y%m%d_%H%M%S)
result_root=${3:-"$repo_root/results/walltime_four_systems_$timestamp"}
timeout_seconds=${COMPARE_TIMEOUT_SECONDS:-1800}
build_jobs=${COMPARE_BUILD_JOBS:-4}
dataset_filter=${COMPARE_DATASETS:-}
fail_on_error=${COMPARE_FAIL_ON_ERROR:-0}

if [[ ! $timeout_seconds =~ ^[1-9][0-9]*$ ]]; then
  echo "COMPARE_TIMEOUT_SECONDS must be a positive integer." >&2
  exit 2
fi
if [[ ! $build_jobs =~ ^[1-9][0-9]*$ ]]; then
  echo "COMPARE_BUILD_JOBS must be a positive integer." >&2
  exit 2
fi
if [[ $fail_on_error != 0 && $fail_on_error != 1 ]]; then
  echo "COMPARE_FAIL_ON_ERROR must be 0 or 1." >&2
  exit 2
fi
if [[ ! -d $adjacency_dir ]]; then
  echo "Adjacency-list directory does not exist: $adjacency_dir" >&2
  exit 2
fi
if [[ ! -d $edge_dir ]]; then
  echo "Edge-list directory does not exist: $edge_dir" >&2
  exit 2
fi
if [[ $(uname -s) != Linux ]]; then
  echo "This experiment requires Linux and GNU coreutils." >&2
  exit 2
fi
for tool in cmake timeout /usr/bin/time awk sha256sum uname find sort date; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool is unavailable: $tool" >&2
    exit 2
  fi
done

adjacency_root=$(cd "$adjacency_dir" && pwd -P)
edge_root=$(cd "$edge_dir" && pwd -P)
mkdir -p "$result_root/build_logs" "$result_root/logs" "$result_root/tomita_inputs"

build_target() {
  local name=$1
  local source_dir=$2
  local build_dir=$3
  local binary_name=$4
  shift 4
  local configure_log="$result_root/build_logs/$name.configure.log"
  local build_log="$result_root/build_logs/$name.build.log"
  echo "CONFIGURE name=$name"
  if ! cmake -S "$source_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release "$@" >"$configure_log" 2>&1; then
    echo "Configure failed; see $configure_log" >&2
    return 1
  fi
  echo "BUILD name=$name"
  if ! cmake --build "$build_dir" --parallel "$build_jobs" >"$build_log" 2>&1; then
    echo "Build failed; see $build_log" >&2
    return 1
  fi
  if [[ ! -x $build_dir/$binary_name ]]; then
    echo "Expected binary was not built: $build_dir/$binary_name" >&2
    return 1
  fi
}

reorder_build="$result_root/build/reorder"
hbbmc_build="$result_root/build/hbbmc"
comparison_build="$result_root/build/comparison_tools"
build_target reorder "$repo_root/our" "$reorder_build" adaptive_bk -DBUILD_TESTING=OFF || exit 2
build_target hbbmc "$repo_root/compare/HBBMCPaperFaithful" "$hbbmc_build" hbbmc_faithful -DBUILD_TESTING=OFF || exit 2
build_target comparison_tools "$repo_root/experiments/walltime" "$comparison_build" tomita_walltime || exit 2

reorder_binary="$reorder_build/adaptive_bk"
hbbmc_binary="$hbbmc_build/hbbmc_faithful"
tomita_adapter="$comparison_build/tomita_stream_adapter"
tomita_binary="$comparison_build/tomita_walltime"
ccrmce_binary="$comparison_build/ccrmce_walltime"
for binary in "$reorder_binary" "$hbbmc_binary" "$tomita_adapter" "$tomita_binary" "$ccrmce_binary"; do
  if [[ ! -x $binary ]]; then
    echo "Missing executable: $binary" >&2
    exit 2
  fi
done

selected_dataset() {
  local group=$1
  local graph=$2
  local relative=$3
  local stem=${graph%.*}
  local wanted
  [[ -z $dataset_filter ]] && return 0
  for wanted in ${dataset_filter//,/ }; do
    if [[ $wanted == "$group" || $wanted == "$graph" ||
          $wanted == "$stem" || $wanted == "$relative" ]]; then
      return 0
    fi
  done
  return 1
}

dataset_specs=()
skipped_graphs=0
while IFS= read -r -d '' adjacency_input; do
  relative=${adjacency_input#"$adjacency_root"/}
  group=${relative%/*}
  [[ $group == "$relative" ]] && group=.
  graph=${relative##*/}
  selected_dataset "$group" "$graph" "$relative" || continue
  if [[ $relative == *","* || $relative == *"|"* ]]; then
    echo "SKIP graph=$relative reason=unsupported_delimiter_in_path" >&2
    skipped_graphs=$((skipped_graphs + 1))
    continue
  fi
  edge_input="$edge_root/$relative"
  if [[ ! -f $edge_input ]]; then
    echo "SKIP graph=$relative reason=missing_edge_pair expected=$edge_input" >&2
    skipped_graphs=$((skipped_graphs + 1))
    continue
  fi
  if ! read -r n m _ <"$adjacency_input" ||
      [[ ! $n =~ ^[0-9]+$ || ! $m =~ ^[0-9]+$ ]]; then
    echo "SKIP graph=$relative reason=invalid_adjacency_header" >&2
    skipped_graphs=$((skipped_graphs + 1))
    continue
  fi
  dataset_specs+=("$relative|$graph|$n|$m|$adjacency_input|$edge_input")
done < <(find "$adjacency_root" -type f ! -path '*/.*' -print0 | sort -z)

if [[ ${#dataset_specs[@]} -eq 0 ]]; then
  echo "No valid paired graph inputs matched." >&2
  exit 2
fi

runs_csv="$result_root/runs.csv"
comparison_csv="$result_root/comparison.csv"
runs_header='variant,relative_path,graph,n,m,status,exit_code,cliques,stored_cliques,wall_time_us,algorithm_time_ms,max_rss_kb'
comparison_header='relative_path,graph,n,m,reorder_status,reorder_cliques,reorder_wall_time_us,hbbmc_status,hbbmc_cliques,hbbmc_wall_time_us,tomita_status,tomita_cliques,tomita_wall_time_us,ccrmce_status,ccrmce_cliques,ccrmce_wall_time_us,all_counts_match'

initialize_csv() {
  local file=$1
  local header=$2
  if [[ ! -f $file ]]; then
    printf '%s\n' "$header" >"$file"
    return
  fi
  local existing
  IFS= read -r existing <"$file" || true
  if [[ $existing != "$header" ]]; then
    echo "Cannot resume $file: CSV header does not match this runner." >&2
    exit 2
  fi
}
initialize_csv "$runs_csv" "$runs_header"
initialize_csv "$comparison_csv" "$comparison_header"

{
  echo "campaign=four-system complete-process wall time"
  echo "repository=$repo_root"
  echo "adjacency_root=$adjacency_root"
  echo "edge_root=$edge_root"
  echo "result_root=$result_root"
  echo "graphs=${#dataset_specs[@]}"
  echo "skipped_graphs=$skipped_graphs"
  echo "minimum_clique_size=3"
  echo "reorder_budget=1000"
  echo "hbbmc_configuration=RMCE graph reduction plus ET3"
  echo "timeout_seconds_per_program=$timeout_seconds"
  echo "wall_time_unit=microseconds"
  echo "wall_time_scope=entire native program invocation including parse, preprocessing, enumeration, retained output work, and teardown"
  echo "tomita_conversion=logged preparation outside the native program wall-time interval"
  echo "execution=sequential"
  sha256sum "$reorder_binary" "$hbbmc_binary" "$tomita_adapter" "$tomita_binary" "$ccrmce_binary"
  uname -a
} >"$result_root/environment.txt"

output_value() {
  local file=$1
  local key=$2
  awk -F= -v key="$key" '$1 == key { value = substr($0, length($1) + 2) } END { print value }' "$file"
}

colon_value() {
  local file=$1
  local key=$2
  awk -F: -v key="$key" '$1 == key { value = substr($0, length($1) + 2) } END { print value }' "$file"
}

already_recorded() {
  local variant=$1
  local relative=$2
  awk -F, -v variant="$variant" -v relative="$relative" 'NR > 1 && $1 == variant && $2 == relative { found = 1 } END { exit !found }' "$runs_csv"
}

comparison_recorded() {
  local relative=$1
  awk -F, -v relative="$relative" 'NR > 1 && $1 == relative { found = 1 } END { exit !found }' "$comparison_csv"
}

write_command() {
  local file=$1
  shift
  printf '%q ' "$@" >"$file"
  printf '\n' >>"$file"
}

campaign_failed=0
run_one() {
  local variant=$1
  local relative=$2
  local graph=$3
  local n=$4
  local m=$5
  shift 5
  local -a command=("$@")
  if already_recorded "$variant" "$relative"; then
    echo "SKIP variant=$variant graph=$relative reason=recorded"
    return
  fi

  local relative_dir safe_graph log_dir base stdout_file stderr_file resource_file
  relative_dir=$(dirname -- "$relative")
  safe_graph=${graph//[^A-Za-z0-9_.-]/_}
  log_dir="$result_root/logs"
  [[ $relative_dir == . ]] || log_dir="$log_dir/$relative_dir"
  mkdir -p "$log_dir"
  base="$log_dir/${safe_graph}__${variant}"
  stdout_file="$base.stdout"
  stderr_file="$base.stderr"
  resource_file="$base.resources"
  write_command "$base.command" "${command[@]}"

  local start_ns stop_ns run_pid heartbeat_pid now_ns exit_code wall_us
  echo "START variant=$variant graph=$relative timeout=${timeout_seconds}s"
  start_ns=$(date +%s%N)
  /usr/bin/time -f 'max_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S\ntime_exit_code=%x' -o "$resource_file" timeout --signal=TERM --kill-after=10s "${timeout_seconds}s" "${command[@]}" >"$stdout_file" 2>"$stderr_file" &
  run_pid=$!
  (
    while true; do
      sleep 30
      kill -0 "$run_pid" 2>/dev/null || exit
      now_ns=$(date +%s%N)
      echo "HEARTBEAT variant=$variant graph=$relative elapsed=$(((now_ns - start_ns) / 1000000000))s"
    done
  ) &
  heartbeat_pid=$!
  wait "$run_pid"
  exit_code=$?
  stop_ns=$(date +%s%N)
  kill "$heartbeat_pid" 2>/dev/null || true
  wait "$heartbeat_pid" 2>/dev/null || true
  wall_us=$(((stop_ns - start_ns) / 1000))

  local cliques= stored= algorithm_time_ms= contract_ok=0
  case "$variant" in
    reorder)
      cliques=$(output_value "$stdout_file" reorder.cliques)
      stored=$(output_value "$stdout_file" reorder.stored_cliques)
      algorithm_time_ms=$(output_value "$stdout_file" reorder.runtime_ms)
      if [[ $stored == "$cliques" &&
            $(output_value "$stdout_file" reorder.minimum_clique_size) == 3 &&
            $(output_value "$stdout_file" reorder.budget) == 1000 ]]; then
        contract_ok=1
      fi
      ;;
    hbbmc)
      cliques=$(output_value "$stdout_file" maximal_cliques)
      algorithm_time_ms=$(output_value "$stdout_file" algorithm_runtime_ms)
      if [[ $(output_value "$stdout_file" algorithm) == HBBMC++ &&
            $(output_value "$stdout_file" graph_reduction) == rmce &&
            $(output_value "$stdout_file" minimum_clique_size) == 3 ]]; then
        contract_ok=1
      fi
      ;;
    tomita)
      cliques=$(output_value "$stdout_file" maximal_cliques)
      stored=$(output_value "$stdout_file" stored_cliques)
      algorithm_time_ms=$(output_value "$stdout_file" algorithm_wall_ms)
      if [[ $stored == "$cliques" &&
            $(output_value "$stdout_file" algorithm) == tomita-adjacency-list &&
            $(output_value "$stdout_file" minimum_clique_size) == 3 &&
            $(output_value "$stdout_file" clique_storage) == all ]]; then
        contract_ok=1
      fi
      ;;
    ccrmce)
      cliques=$(colon_value "$stdout_file" Mclique)
      stored=$(colon_value "$stdout_file" stored_cliques)
      algorithm_time_ms=$(colon_value "$stdout_file" time)
      algorithm_time_ms=${algorithm_time_ms%ms}
      if [[ $stored == "$cliques" &&
            $(colon_value "$stdout_file" clique_storage) == all ]]; then
        contract_ok=1
      fi
      ;;
  esac

  local status max_rss
  max_rss=$(output_value "$resource_file" max_rss_kb)
  if [[ $exit_code -eq 0 && $contract_ok -eq 1 &&
        $cliques =~ ^[0-9]+$ &&
        $algorithm_time_ms =~ ^[0-9]+([.][0-9]+)?$ &&
        $max_rss =~ ^[0-9]+$ ]]; then
    status=completed
  elif [[ $exit_code -eq 124 || $exit_code -eq 137 ]]; then
    status=timeout
    campaign_failed=1
  elif [[ $exit_code -eq 0 ]]; then
    status=invalid_output
    campaign_failed=1
  else
    status=failed
    campaign_failed=1
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$variant" "$relative" "$graph" "$n" "$m" "$status" "$exit_code" "$cliques" "$stored" "$wall_us" "$algorithm_time_ms" "$max_rss" >>"$runs_csv"
  echo "DONE variant=$variant graph=$relative status=$status wall_time_us=$wall_us cliques=${cliques:-NA}"
}

record_tomita_input_failure() {
  local relative=$1
  local graph=$2
  local n=$3
  local m=$4
  already_recorded tomita "$relative" && return
  printf 'tomita,%s,%s,%s,%s,input_error,1,,,,,\n' "$relative" "$graph" "$n" "$m" >>"$runs_csv"
  campaign_failed=1
}

csv_field() {
  local variant=$1
  local relative=$2
  local field=$3
  awk -F, -v variant="$variant" -v relative="$relative" -v field="$field" '$1 == variant && $2 == relative { value = $field } END { print value }' "$runs_csv"
}

append_comparison() {
  local relative=$1
  local graph=$2
  local n=$3
  local m=$4
  comparison_recorded "$relative" && return
  local r_status r_count r_wall h_status h_count h_wall
  local t_status t_count t_wall c_status c_count c_wall all_counts=NA
  r_status=$(csv_field reorder "$relative" 6)
  r_count=$(csv_field reorder "$relative" 8)
  r_wall=$(csv_field reorder "$relative" 10)
  h_status=$(csv_field hbbmc "$relative" 6)
  h_count=$(csv_field hbbmc "$relative" 8)
  h_wall=$(csv_field hbbmc "$relative" 10)
  t_status=$(csv_field tomita "$relative" 6)
  t_count=$(csv_field tomita "$relative" 8)
  t_wall=$(csv_field tomita "$relative" 10)
  c_status=$(csv_field ccrmce "$relative" 6)
  c_count=$(csv_field ccrmce "$relative" 8)
  c_wall=$(csv_field ccrmce "$relative" 10)
  if [[ $r_status == completed && $h_status == completed &&
        $t_status == completed && $c_status == completed ]]; then
    if [[ $r_count == "$h_count" && $r_count == "$t_count" &&
          $r_count == "$c_count" ]]; then
      all_counts=yes
    else
      all_counts=no
      campaign_failed=1
    fi
  else
    campaign_failed=1
  fi
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$relative" "$graph" "$n" "$m" "$r_status" "$r_count" "$r_wall" "$h_status" "$h_count" "$h_wall" "$t_status" "$t_count" "$t_wall" "$c_status" "$c_count" "$c_wall" "$all_counts" >>"$comparison_csv"
  echo "COMPARE graph=$relative all_counts_match=$all_counts"
}

for spec in "${dataset_specs[@]}"; do
  IFS='|' read -r relative graph n m adjacency_input edge_input <<<"$spec"
  relative_dir=$(dirname -- "$relative")
  tomita_dir="$result_root/tomita_inputs"
  [[ $relative_dir == . ]] || tomita_dir="$tomita_dir/$relative_dir"
  mkdir -p "$tomita_dir"
  tomita_input="$tomita_dir/${graph}.input"
  adapter_base="$tomita_dir/${graph}.adapter"

  if ! already_recorded tomita "$relative"; then
    write_command "$adapter_base.command" "$tomita_adapter" "$edge_input" "$n" "$m"
    if ! "$tomita_adapter" "$edge_input" "$n" "$m" >"$tomita_input" 2>"$adapter_base.stderr"; then
      echo "Tomita input conversion failed for $relative; recording the failure and continuing." >&2
      record_tomita_input_failure "$relative" "$graph" "$n" "$m"
    fi
  fi

  run_one reorder "$relative" "$graph" "$n" "$m" env OMP_NUM_THREADS=1 "$reorder_binary" "$adjacency_input" --budget 1000 --min-clique-size 3
  run_one hbbmc "$relative" "$graph" "$n" "$m" env OMP_NUM_THREADS=1 "$hbbmc_binary" "$edge_input" --graph-reduction rmce --et 3 --num-vertices "$n" --min-clique-size 3
  run_one tomita "$relative" "$graph" "$n" "$m" env OMP_NUM_THREADS=1 "$tomita_binary" "$tomita_input"
  run_one ccrmce "$relative" "$graph" "$n" "$m" env OMP_NUM_THREADS=1 "$ccrmce_binary" noUVM -f_txt "$edge_input"
  append_comparison "$relative" "$graph" "$n" "$m"
done

if ! awk -F, 'NR > 1 && $6 != "completed" { exit 1 }' "$runs_csv"; then
  campaign_failed=1
fi
if ! awk -F, 'NR > 1 && $17 != "yes" { exit 1 }' "$comparison_csv"; then
  campaign_failed=1
fi

echo "CAMPAIGN_COMPLETE graphs=${#dataset_specs[@]} skipped=$skipped_graphs results=$result_root"
if [[ $fail_on_error == 1 && $campaign_failed -ne 0 ]]; then
  exit 1
fi
