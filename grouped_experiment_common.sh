#!/usr/bin/env bash

# Shared discovery helpers for the grouped graph corpus. Callers source
# experiment_common.sh first so graph-header and validation helpers are present.

final_graph_filter_matches() {
  local group=$1
  local graph_file=$2
  local relative_path=$3
  local filter=$4
  local graph_stem=${graph_file%.*}
  local wanted wanted_key stem_key

  [[ -z $filter ]] && return 0
  stem_key=${graph_stem,,}
  stem_key=${stem_key//[^[:alnum:]]/}
  for wanted in ${filter//,/ }; do
    if [[ $wanted == "$group" || $wanted == "$graph_file" ||
          $wanted == "$graph_stem" || $wanted == "$relative_path" ]]; then
      return 0
    fi
    wanted_key=${wanted,,}
    wanted_key=${wanted_key//[^[:alnum:]]/}
    [[ -n $wanted_key && $wanted_key == "$stem_key" ]] && return 0
  done
  return 1
}

# Discover ADJACENCY_ROOT/GROUP/GRAPH inputs. A paired edge-list input is the
# file at the identical relative path below EDGE_ROOT. Direct files are
# supported as group "." for small fixtures.
#
# Populates FINAL_SELECTED_GROUPS and FINAL_SELECTED_DATASETS. Each dataset
# entry is GROUP|GRAPH_FILE|N|M|PURE_INPUT|EDGE_INPUT. A malformed adjacency
# file or missing pair is reported and skipped so one bad graph cannot abort a
# campaign that still has runnable inputs.
final_collect_grouped_datasets() {
  local adjacency_root=$1
  local edge_root=$2
  local filter=$3
  local require_edge=$4
  local pure_input relative_path group graph_file hbbmc_input edge_input
  declare -A seen_groups=()

  if [[ ! -d $adjacency_root ]]; then
    echo "Adjacency-list directory does not exist: $adjacency_root" >&2
    return 1
  fi
  FINAL_ADJACENCY_ROOT=$(cd "$adjacency_root" && pwd -P)
  FINAL_EDGE_ROOT=
  if [[ $require_edge == 1 && ! -d $edge_root ]]; then
    echo "Edge-list directory does not exist: $edge_root" >&2
    return 1
  fi
  if [[ -n $edge_root ]]; then
    FINAL_EDGE_ROOT=$(cd "$edge_root" && pwd -P)
  fi

  FINAL_SELECTED_GROUPS=()
  FINAL_SELECTED_DATASETS=()
  FINAL_SKIPPED_DATASETS=0
  while IFS= read -r -d '' pure_input; do
    relative_path=${pure_input#"$FINAL_ADJACENCY_ROOT"/}
    group=${relative_path%/*}
    [[ $group == "$relative_path" ]] && group=.
    graph_file=${relative_path##*/}
    final_graph_filter_matches "$group" "$graph_file" "$relative_path" \
      "$filter" || continue

    hbbmc_input=
    if [[ -n $FINAL_EDGE_ROOT ]]; then
      hbbmc_input="$FINAL_EDGE_ROOT/$relative_path"
    fi
    if [[ $require_edge == 1 && ! -f $hbbmc_input ]]; then
      echo "SKIP graph=$relative_path reason=missing_edge_pair expected=$hbbmc_input" >&2
      FINAL_SKIPPED_DATASETS=$((FINAL_SKIPPED_DATASETS + 1))
      continue
    fi
    if ! final_read_graph_header "$pure_input"; then
      echo "SKIP graph=$relative_path reason=invalid_adjacency_header" >&2
      FINAL_SKIPPED_DATASETS=$((FINAL_SKIPPED_DATASETS + 1))
      continue
    fi
    if [[ -z ${seen_groups[$group]:-} ]]; then
      FINAL_SELECTED_GROUPS+=("$group")
      seen_groups[$group]=1
    fi
    FINAL_SELECTED_DATASETS+=(
      "$group|$graph_file|$FINAL_GRAPH_N|$FINAL_GRAPH_M|$pure_input|$hbbmc_input"
    )
  done < <(find "$FINAL_ADJACENCY_ROOT" -type f ! -path '*/.*' -print0 | sort -z)

  if [[ ${#FINAL_SELECTED_DATASETS[@]} -eq 0 ]]; then
    echo "No graph files matched under $FINAL_ADJACENCY_ROOT." >&2
    return 1
  fi

}

final_group_result_dir() {
  local result_root=$1
  local group=$2
  if [[ $group == . ]]; then
    printf '%s' "$result_root"
  else
    printf '%s/%s' "$result_root" "$group"
  fi
}
