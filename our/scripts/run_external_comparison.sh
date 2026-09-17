#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage:
  run_external_comparison.sh DATA_DIR OUTPUT_CSV TIMEOUT_SECONDS REPETITIONS \
      BUDGET MIN_CLIQUE_SIZE DATASETS LABEL=BINARY [LABEL=BINARY ...]

DATASETS is "all" or a comma-separated list of basenames, with or without
".txt". Existing OUTPUT_CSV files are never overwritten.

Example:
  ./scripts/run_external_comparison.sh /graphs results/runs.csv 30 1 1000 3 \
      all ccrmce=./build/adaptive_bk pxr_no_et=/tmp/pxr/adaptive_bk
USAGE
}

if (( $# < 8 )); then
  usage
  exit 2
fi

data_dir=$1
output_csv=$2
timeout_seconds=$3
repetitions=$4
budget=$5
minimum_size=$6
dataset_spec=$7
shift 7
variant_specs=("$@")

if [[ ! -d $data_dir ]]; then
  echo "dataset directory does not exist: $data_dir" >&2
  exit 2
fi
if [[ -e $output_csv ]]; then
  echo "refusing to overwrite existing output: $output_csv" >&2
  exit 2
fi
if [[ ! $timeout_seconds =~ ^[1-9][0-9]*$ ||
      ! $repetitions =~ ^[1-9][0-9]*$ ||
      ! $budget =~ ^[0-9]+$ ||
      ! $minimum_size =~ ^[1-9][0-9]*$ ]]; then
  echo "timeout/repetitions/minimum size must be positive integers; budget must be nonnegative" >&2
  exit 2
fi

variant_labels=()
variant_binaries=()
for specification in "${variant_specs[@]}"; do
  if [[ $specification != *=* ]]; then
    echo "variant must have LABEL=BINARY form: $specification" >&2
    exit 2
  fi
  label=${specification%%=*}
  binary=${specification#*=}
  if [[ ! $label =~ ^[A-Za-z0-9_.-]+$ ]]; then
    echo "invalid variant label: $label" >&2
    exit 2
  fi
  if [[ ! -x $binary ]]; then
    echo "variant binary is not executable: $binary" >&2
    exit 2
  fi
  variant_labels+=("$label")
  variant_binaries+=("$binary")
done

graphs=()
if [[ $dataset_spec == all ]]; then
  shopt -s nullglob
  graphs=("$data_dir"/*.txt)
  shopt -u nullglob
else
  IFS=',' read -r -a requested_datasets <<< "$dataset_spec"
  for name in "${requested_datasets[@]}"; do
    [[ $name == *.txt ]] || name+=.txt
    graphs+=("$data_dir/$name")
  done
fi
if (( ${#graphs[@]} == 0 )); then
  echo "no datasets selected" >&2
  exit 2
fi
for graph in "${graphs[@]}"; do
  if [[ ! -f $graph ]]; then
    echo "dataset does not exist: $graph" >&2
    exit 2
  fi
done

output_dir=$(dirname -- "$output_csv")
mkdir -p -- "$output_dir"
scratch_dir=$(mktemp -d /tmp/ccrmce_benchmark.XXXXXX)
cleanup() {
  rm -rf -- "$scratch_dir"
}
trap cleanup EXIT

printf '%s\n' 'dataset,variant,repetition,budget,min_clique_size,status,exit_code,cliques,runtime_ms,wall_seconds,max_rss_kb' > "$output_csv"

declare -A expected_counts
completed=0
timed_out=0
failed=0
mismatched=0

for graph in "${graphs[@]}"; do
  dataset=$(basename -- "$graph")
  for (( repetition = 1; repetition <= repetitions; ++repetition )); do
    key="$dataset:$repetition"
    for (( variant = 0; variant < ${#variant_labels[@]}; ++variant )); do
      label=${variant_labels[variant]}
      binary=${variant_binaries[variant]}
      stdout_file="$scratch_dir/stdout"
      stderr_file="$scratch_dir/stderr"
      timing_file="$scratch_dir/timing"

      echo "RUN dataset=$dataset variant=$label repetition=$repetition" >&2
      /usr/bin/time -f '%e,%M' -o "$timing_file" \
        timeout --foreground --signal=TERM --kill-after=5 \
          "${timeout_seconds}s" "$binary" "$graph" \
          --budget "$budget" --min-clique-size "$minimum_size" \
          > "$stdout_file" 2> "$stderr_file"
      exit_code=$?

      cliques=$(awk -F= '$1 == "reorder.cliques" { value=$2 } END { print value }' "$stdout_file")
      runtime_ms=$(awk -F= '$1 == "reorder.runtime_ms" { value=$2 } END { print value }' "$stdout_file")
      timing=$(awk -F, 'NF == 2 { value=$0 } END { print value }' "$timing_file")
      wall_seconds=${timing%%,*}
      max_rss_kb=${timing#*,}
      [[ $timing == *,* ]] || {
        wall_seconds=
        max_rss_kb=
      }

      if (( exit_code == 0 )) && [[ $cliques =~ ^[0-9]+$ &&
                                    $runtime_ms =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        status=ok
        if [[ -z ${expected_counts[$key]+present} ]]; then
          expected_counts[$key]=$cliques
        elif [[ ${expected_counts[$key]} != "$cliques" ]]; then
          status=count_mismatch
          (( ++mismatched ))
        fi
        (( ++completed ))
      elif (( exit_code == 124 )); then
        status=timeout
        (( ++timed_out ))
      else
        status=error
        (( ++failed ))
      fi

      printf '%s,%s,%d,%s,%s,%s,%d,%s,%s,%s,%s\n' \
        "$dataset" "$label" "$repetition" "$budget" "$minimum_size" \
        "$status" "$exit_code" "$cliques" "$runtime_ms" \
        "$wall_seconds" "$max_rss_kb" >> "$output_csv"
      echo "DONE dataset=$dataset variant=$label status=$status runtime_ms=${runtime_ms:-NA}" >&2
    done
  done
done

echo "BENCHMARK_COMPLETE completed=$completed timeouts=$timed_out errors=$failed count_mismatches=$mismatched output=$output_csv"
if (( failed != 0 || mismatched != 0 )); then
  exit 1
fi
