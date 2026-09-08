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

# Discover adjacencylist/GROUP/GRAPH inputs. A paired HBBMC input is the file
# at the identical relative path below edgelist/. Direct files are supported
# as group "." for small fixtures.
#
# Populates FINAL_SELECTED_GROUPS and FINAL_SELECTED_DATASETS. Each dataset
# entry is GROUP|GRAPH_FILE|N|M|PURE_INPUT|HBBMC_INPUT.
final_collect_grouped_datasets() {
  local data_root=$1
  local filter=$2
  local require_hbbmc=$3
  local pure_input relative_path group graph_file hbbmc_input edge_input
  declare -A seen_groups=()

  FINAL_ADJACENCY_ROOT="$data_root/adjacencylist"
  FINAL_EDGE_ROOT="$data_root/edgelist"
  if [[ ! -d $FINAL_ADJACENCY_ROOT ]]; then
    echo "Expected an adjacencylist/ directory under $data_root." >&2
    return 1
  fi
  if [[ $require_hbbmc == 1 && ! -d $FINAL_EDGE_ROOT ]]; then
    echo "Expected an edgelist/ directory under $data_root." >&2
    return 1
  fi

  FINAL_SELECTED_GROUPS=()
  FINAL_SELECTED_DATASETS=()
  while IFS= read -r -d '' pure_input; do
    relative_path=${pure_input#"$FINAL_ADJACENCY_ROOT"/}
    group=${relative_path%/*}
    [[ $group == "$relative_path" ]] && group=.
    graph_file=${relative_path##*/}
    final_graph_filter_matches "$group" "$graph_file" "$relative_path" \
      "$filter" || continue

    hbbmc_input="$FINAL_EDGE_ROOT/$relative_path"
    if [[ $require_hbbmc == 1 && ! -f $hbbmc_input ]]; then
      echo "Missing paired HBBMC input: $hbbmc_input" >&2
      return 1
    fi
    final_read_graph_header "$pure_input" || return 1
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

  if [[ $require_hbbmc == 1 ]]; then
    while IFS= read -r -d '' edge_input; do
      relative_path=${edge_input#"$FINAL_EDGE_ROOT"/}
      if [[ ! -f $FINAL_ADJACENCY_ROOT/$relative_path ]]; then
        echo "Missing paired Pure input: $FINAL_ADJACENCY_ROOT/$relative_path" >&2
        return 1
      fi
    done < <(find "$FINAL_EDGE_ROOT" -type f ! -path '*/.*' -print0 | sort -z)
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
