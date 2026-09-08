#!/usr/bin/env bash

# Shared experiment definitions and shell helpers. Every runner supplies its
# own physical source directory and build directory so CMake options cannot
# leak between ablation variants.

FINAL_PURE_BINARY_NAME=pure_mce
FINAL_HBBMC_BINARY_NAME=hbbmc_faithful

FINAL_PRUNING_RULES=(
  normalization
  subsumption
  unit
  usefulness
  antichain
  fail_first
  zero_coverage
)

final_require_positive_integer() {
  local name=$1
  local value=$2
  if [[ ! $value =~ ^[1-9][0-9]*$ ]]; then
    echo "$name must be a positive integer (received: $value)." >&2
    return 1
  fi
}

final_require_nonnegative_integer() {
  local name=$1
  local value=$2
  if [[ ! $value =~ ^[0-9]+$ ]]; then
    echo "$name must be a non-negative integer (received: $value)." >&2
    return 1
  fi
}

final_require_tools() {
  local tool
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || {
      echo "Required tool is unavailable: $tool" >&2
      return 1
    }
  done
}

final_read_graph_header() {
  local input=$1
  if ! read -r FINAL_GRAPH_N FINAL_GRAPH_M _ <"$input" ||
      [[ ! $FINAL_GRAPH_N =~ ^[0-9]+$ ||
         ! $FINAL_GRAPH_M =~ ^[0-9]+$ ]]; then
    echo "Invalid Pure graph header: $input" >&2
    return 1
  fi
}

final_build() {
  local name=$1
  local source_dir=$2
  local build_dir=$3
  local log_root=$4
  local build_jobs=$5
  shift 5
  local -a options=("$@")

  mkdir -p "$build_dir" "$log_root"
  echo "CONFIGURE name=$name build=$build_dir"
  cmake -S "$source_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release "${options[@]}" >"$log_root/${name}.configure.log" 2>&1 || {
    echo "Configure failed; see $log_root/${name}.configure.log" >&2
    return 1
  }
  echo "BUILD name=$name"
  cmake --build "$build_dir" --parallel "$build_jobs" \
    >"$log_root/${name}.build.log" 2>&1 || {
    echo "Build failed; see $log_root/${name}.build.log" >&2
    return 1
  }
}

final_write_command() {
  local file=$1
  shift
  printf '%q ' "$@" >"$file"
  printf '\n' >>"$file"
}

final_output_value() {
  local file=$1
  local key=$2
  awk -F= -v key="$key" '$1 == key {
      value = substr($0, length($1) + 2)
    }
    END { print value }' "$file"
}

final_pruning_config_matches() {
  local output=$1
  local disabled=${2:-none}
  local rule expected actual
  for rule in "${FINAL_PRUNING_RULES[@]}"; do
    expected=1
    [[ $rule == "$disabled" ]] && expected=0
    actual=$(final_output_value "$output" "pure.config.pruning.$rule")
    [[ $actual == "$expected" ]] || return 1
  done
}

final_all_uint() {
  local value
  for value in "$@"; do
    [[ $value =~ ^[0-9]+$ ]] || return 1
  done
}

final_all_number() {
  local value
  for value in "$@"; do
    [[ $value =~ ^[0-9]+([.][0-9]+)?$ ]] || return 1
  done
}

final_run_timed() {
  local label=$1
  local stdout_file=$2
  local stderr_file=$3
  local resource_file=$4
  local timeout_seconds=$5
  shift 5
  local start_ns stop_ns run_pid heartbeat_pid now_ns

  echo "START $label timeout=${timeout_seconds}s"
  start_ns=$(date +%s%N)
  /usr/bin/time \
    -f 'max_rss_kb=%M\nuser_seconds=%U\nsystem_seconds=%S\ntime_exit_code=%x' \
    -o "$resource_file" \
    timeout --signal=TERM --kill-after=10s "${timeout_seconds}s" \
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
  FINAL_RUN_EXIT_CODE=$?
  stop_ns=$(date +%s%N)
  kill "$heartbeat_pid" 2>/dev/null || true
  wait "$heartbeat_pid" 2>/dev/null || true
  FINAL_RUN_WALL_MS=$(((stop_ns - start_ns) / 1000000))
}

final_status_from_exit() {
  local exit_code=$1
  if [[ $exit_code -eq 124 || $exit_code -eq 137 ]]; then
    printf 'timeout'
  elif [[ $exit_code -eq 0 ]]; then
    printf 'invalid_output'
  else
    printf 'failed'
  fi
}
