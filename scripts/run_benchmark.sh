#!/usr/bin/env bash

set -uo pipefail

usage() {
  cat <<'EOF'
Usage: ./run_compare [DATA_ROOT [RESULT_ROOT]]

Build the current repository and compare:
  - independent paper-faithful HBBMC++ (RMCE + ET(3)); and
  - Pure ReorderSib mode 1 with its remaining defaults.

Defaults:
  DATA_ROOT    /data/labdata/wajid/hbbmcData
  RESULT_ROOT  REPO_ROOT/results/reorder_mode1_vs_hbbmc_faithful_TIMESTAMP
  timeout      600 seconds per algorithm/dataset run

Builds (Release, stripped, in-tree — Linux/Ubuntu toolchain required):
  our/build/bk_algorithm
  compare/HBBMCPaperFaithful/build/hbbmc_faithful

Environment overrides:
  COMPARE_TIMEOUT_SECONDS  timeout per run (default: 600)
  COMPARE_BUILD_JOBS       parallel build jobs (default: all cores)
  COMPARE_DATASETS         comma/space-separated dataset names to run

The result directory contains runs.csv, comparison.csv, environment.txt,
build logs, and stdout/stderr/resource logs for every run. Existing rows
are skipped when an explicit RESULT_ROOT is resumed.
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" || $# -gt 2 ]]; then
  usage
  [[ $# -gt 2 ]] && exit 2
  exit 0
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
data_root=${1:-/data/labdata/wajid/hbbmcData}
timestamp=$(date +%Y%m%d_%H%M%S)
result_root=${2:-"$repo_root/results/reorder_mode1_vs_hbbmc_faithful_$timestamp"}
timeout_seconds=${COMPARE_TIMEOUT_SECONDS:-600}
build_jobs=${COMPARE_BUILD_JOBS:-}
dataset_filter=${COMPARE_DATASETS:-}

if [[ $(uname -s) != Linux ]]; then
  echo "This build targets Ubuntu/Linux (uses GNU coreutils: sha256sum, GNU time, GNU timeout, GNU strip)." >&2
  echo "Detected: $(uname -s). Run this on the Ubuntu host." >&2
  exit 2
fi

if [[ ! $timeout_seconds =~ ^[1-9][0-9]*$ ]]; then
  echo "COMPARE_TIMEOUT_SECONDS must be a positive integer." >&2
  exit 2
fi
if [[ -n $build_jobs && ! $build_jobs =~ ^[1-9][0-9]*$ ]]; then
  echo "COMPARE_BUILD_JOBS must be a positive integer." >&2
  exit 2
fi
if [[ ! -d "$data_root/hbbmc" || ! -d "$data_root/pure" ]]; then
  echo "Expected paired input directories at:" >&2
  echo "  $data_root/hbbmc" >&2
  echo "  $data_root/pure" >&2
  exit 2
fi
for tool in cmake timeout /usr/bin/time awk sed sha256sum strip; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool is unavailable: $tool" >&2
    exit 2
  fi
done

mkdir -p "$result_root/build_logs" "$result_root/logs"

build_target() {
  local name=$1
  local source_dir=$2
  local binary_name=$3
  local target_build_dir="$source_dir/build"
  local configure_log="$result_root/build_logs/${name}_configure.log"
  local build_log="$result_root/build_logs/${name}_build.log"
  local -a build_command=(cmake --build "$target_build_dir")

  if [[ -n $build_jobs ]]; then
    build_command+=(--parallel "$build_jobs")
  else
    build_command+=(--parallel)
  fi

  echo "CONFIGURE name=$name source=$source_dir build=$target_build_dir"
  if ! cmake -S "$source_dir" -B "$target_build_dir" \
      -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
      -DREORDERSIB_PROFILING=OFF -DREORDERSIB_SOLVER_TRACE=OFF \
      >"$configure_log" 2>&1; then
    echo "CMake configure failed; see $configure_log" >&2
    return 1
  fi

  echo "BUILD name=$name"
  if ! "${build_command[@]}" >"$build_log" 2>&1; then
    echo "Build failed; see $build_log" >&2
    return 1
  fi
  if [[ ! -x "$target_build_dir/$binary_name" ]]; then
    echo "Expected binary was not built: $target_build_dir/$binary_name" >&2
    return 1
  fi

  strip --strip-unneeded "$target_build_dir/$binary_name"
}

build_target pure "$repo_root/our" bk_algorithm || exit 2
build_target hbbmc_faithful "$repo_root/compare/HBBMCPaperFaithful" \
  hbbmc_faithful || exit 2

pure_binary="$repo_root/our/build/bk_algorithm"
hbbmc_binary="$repo_root/compare/HBBMCPaperFaithful/build/hbbmc_faithful"
runs_csv="$result_root/runs.csv"
comparison_csv="$result_root/comparison.csv"

if [[ ! -f $runs_csv ]]; then
  printf '%s\n' \
    'variant,dataset,n,m,status,exit_code,total_clique_count,wall_time_ms,wall_time_seconds,algorithm_time_ms' \
    >"$runs_csv"
fi
if [[ ! -f $comparison_csv ]]; then
  printf '%s\n' \
    'dataset,n,m,hbbmc_status,hbbmc_total_cliques,hbbmc_wall_time_ms,hbbmc_algorithm_time_ms,reorder_status,reorder_total_cliques,reorder_wall_time_ms,reorder_algorithm_time_ms,count_match' \
    >"$comparison_csv"
fi

{
  echo "campaign=paper-faithful HBBMC++ versus Pure ReorderSib mode 1"
  echo "repository=$repo_root"
  echo "data_root=$data_root"
  echo "result_root=$result_root"
  echo "minimum_clique_size=3"
  echo "timeout_seconds_per_run=$timeout_seconds"
  echo "execution=sequential"
  echo "pure_command=bk_algorithm GRAPH 1"
  echo "pure_defaults=order1,minCliqueSize3,budget10000"
  echo "hbbmc_command=hbbmc_faithful GRAPH --graph-reduction rmce --et 3 --num-vertices N --min-clique-size 3"
  sha256sum "$pure_binary" "$hbbmc_binary"
  uname -a
} >"$result_root/environment.txt"

already_recorded() {
  local variant=$1
  local dataset=$2
  awk -F, -v variant="$variant" -v dataset="$dataset" \
    'NR > 1 && $1 == variant && $2 == dataset { found = 1 } END { exit !found }' \
    "$runs_csv"
}

comparison_recorded() {
  local dataset=$1
  awk -F, -v dataset="$dataset" \
    'NR > 1 && $1 == dataset { found = 1 } END { exit !found }' \
    "$comparison_csv"
}

selected_dataset() {
  local dataset=$1
  local hbbmc_name=$2
  local pure_name=$3
  local wanted

  [[ -z $dataset_filter ]] && return 0
  for wanted in ${dataset_filter//,/ }; do
    if [[ $wanted == "$dataset" || $wanted == "$hbbmc_name" ||
          $wanted == "$pure_name" ]]; then
      return 0
    fi
  done
  return 1
}

csv_field() {
  local variant=$1
  local dataset=$2
  local field=$3
  awk -F, -v variant="$variant" -v dataset="$dataset" -v field="$field" \
    '$1 == variant && $2 == dataset { value = $field } END { print value }' \
    "$runs_csv"
}

run_one() {
  local variant=$1
  local dataset=$2
  local n=$3
  local m=$4
  local input_path=$5
  shift 5
  local -a command=("$@")

  if already_recorded "$variant" "$dataset"; then
    echo "SKIP variant=$variant dataset=$dataset reason=already_recorded"
    return
  fi

  local safe_dataset=${dataset//[^A-Za-z0-9_.-]/_}
  local log_base="$result_root/logs/${variant}__${safe_dataset}"
  local stdout_file="$log_base.stdout"
  local stderr_file="$log_base.stderr"
  local resource_file="$log_base.resources"
  local command_file="$log_base.command"
  local start_ns stop_ns run_pid heartbeat_pid exit_code wall_time_ms
  local wall_time_seconds total_cliques algorithm_time_ms status

  {
    printf 'input=%q ' "$input_path"
    printf '%q ' "${command[@]}"
    printf '\n'
  } >"$command_file"

  echo "START variant=$variant dataset=$dataset timeout=${timeout_seconds}s"
  start_ns=$(date +%s%N)
  /usr/bin/time -f 'max_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S\ntime_exit_code=%x' \
    -o "$resource_file" \
    timeout --signal=TERM --kill-after=10s "${timeout_seconds}s" \
    "${command[@]}" >"$stdout_file" 2>"$stderr_file" &
  run_pid=$!
  (
    while true; do
      sleep 30
      if ! kill -0 "$run_pid" 2>/dev/null; then
        exit
      fi
      local_now_ns=$(date +%s%N)
      echo "HEARTBEAT variant=$variant dataset=$dataset elapsed=$(((local_now_ns - start_ns) / 1000000000))s"
    done
  ) &
  heartbeat_pid=$!

  wait "$run_pid"
  exit_code=$?
  stop_ns=$(date +%s%N)
  kill "$heartbeat_pid" 2>/dev/null || true
  wait "$heartbeat_pid" 2>/dev/null || true
  wall_time_ms=$(((stop_ns - start_ns) / 1000000))
  wall_time_seconds=$(awk -v ms="$wall_time_ms" \
    'BEGIN { printf "%.3f", ms / 1000.0 }')

  total_cliques=""
  algorithm_time_ms=""
  if [[ $variant == reorder_mode1 ]]; then
    local summary
    summary=$(sed -n '/^PureReorderSib:/p' "$stdout_file" | tail -1)
    total_cliques=$(sed -n \
      's/.*cliques=\([0-9][0-9]*\).*/\1/p' <<<"$summary")
    algorithm_time_ms=$(sed -n \
      's/.*time=\([0-9][0-9.]*\) ms.*/\1/p' <<<"$summary")
  else
    total_cliques=$(sed -n 's/^maximal_cliques=//p' "$stdout_file" | tail -1)
    algorithm_time_ms=$(sed -n 's/^algorithm_runtime_ms=//p' "$stdout_file" | tail -1)
  fi

  if [[ $exit_code -eq 0 && -n $total_cliques && -n $algorithm_time_ms ]]; then
    status=completed
  elif [[ $exit_code -eq 124 || $exit_code -eq 137 ]]; then
    status=timeout
  else
    status=failed
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$variant" "$dataset" "$n" "$m" "$status" "$exit_code" \
    "$total_cliques" "$wall_time_ms" "$wall_time_seconds" \
    "$algorithm_time_ms" >>"$runs_csv"

  echo "DONE variant=$variant dataset=$dataset status=$status wall=${wall_time_seconds}s cliques=${total_cliques:-NA} algorithm_ms=${algorithm_time_ms:-NA}"
}

append_comparison() {
  local dataset=$1
  local n=$2
  local m=$3
  local h_status h_count h_wall h_algorithm
  local r_status r_count r_wall r_algorithm count_match

  comparison_recorded "$dataset" && return
  h_status=$(csv_field hbbmc_faithful "$dataset" 5)
  h_count=$(csv_field hbbmc_faithful "$dataset" 7)
  h_wall=$(csv_field hbbmc_faithful "$dataset" 8)
  h_algorithm=$(csv_field hbbmc_faithful "$dataset" 10)
  r_status=$(csv_field reorder_mode1 "$dataset" 5)
  r_count=$(csv_field reorder_mode1 "$dataset" 7)
  r_wall=$(csv_field reorder_mode1 "$dataset" 8)
  r_algorithm=$(csv_field reorder_mode1 "$dataset" 10)

  count_match=NA
  if [[ $h_status == completed && $r_status == completed &&
        -n $h_count && -n $r_count ]]; then
    if [[ $h_count == "$r_count" ]]; then
      count_match=yes
    else
      count_match=no
    fi
  fi

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$dataset" "$n" "$m" "$h_status" "$h_count" "$h_wall" \
    "$h_algorithm" "$r_status" "$r_count" "$r_wall" "$r_algorithm" \
    "$count_match" >>"$comparison_csv"
  echo "COMPARE dataset=$dataset count_match=$count_match"
}

pair_count=0
missing_count=0
while IFS= read -r hbbmc_input; do
  hbbmc_name=$(basename "$hbbmc_input")
  case $hbbmc_name in
    *.edges)
      dataset=${hbbmc_name%.edges}
      pure_name="$dataset.graph"
      ;;
    *.clean)
      pure_name=${hbbmc_name%.clean}
      dataset=${pure_name%.txt}
      ;;
    *)
      echo "SKIP input=$hbbmc_input reason=unrecognized_hbbmc_suffix" >&2
      continue
      ;;
  esac

  if ! selected_dataset "$dataset" "$hbbmc_name" "$pure_name"; then
    continue
  fi

  pure_input="$data_root/pure/$pure_name"
  if [[ ! -f $pure_input ]]; then
    echo "SKIP input=$hbbmc_input reason=missing_pair expected=$pure_input" >&2
    missing_count=$((missing_count + 1))
    continue
  fi
  if ! read -r n m _ <"$pure_input" ||
      [[ ! $n =~ ^[0-9]+$ || ! $m =~ ^[0-9]+$ ]]; then
    echo "SKIP input=$pure_input reason=invalid_pure_header" >&2
    missing_count=$((missing_count + 1))
    continue
  fi

  pair_count=$((pair_count + 1))
  run_one hbbmc_faithful "$dataset" "$n" "$m" "$hbbmc_input" \
    env OMP_NUM_THREADS=1 \
    "$hbbmc_binary" "$hbbmc_input" \
    --graph-reduction rmce --et 3 --num-vertices "$n" --min-clique-size 3

  run_one reorder_mode1 "$dataset" "$n" "$m" "$pure_input" \
    env -u PURE_HITSET_BUDGET -u VLDB_VALIDATION -u VLDB_PRINT_CLIQUES \
    OMP_NUM_THREADS=1 \
    "$pure_binary" "$pure_input" 1

  append_comparison "$dataset" "$n" "$m"
done < <(find "$data_root/hbbmc" -maxdepth 1 -type f -print | LC_ALL=C sort)

if [[ $pair_count -eq 0 ]]; then
  echo "No paired datasets matched${dataset_filter:+ filter '$dataset_filter'}." >&2
  exit 2
fi

echo "CAMPAIGN_COMPLETE pairs=$pair_count missing=$missing_count"
echo "RUNS=$runs_csv"
echo "COMPARISON=$comparison_csv"
