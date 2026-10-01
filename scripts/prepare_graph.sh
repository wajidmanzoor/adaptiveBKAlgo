#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/prepare_graph.sh RAW_EDGE_LIST ADJACENCY_OUTPUT EDGE_OUTPUT [NUM_VERTICES]

Convert a whitespace-separated edge list into the two paired formats used by
the artifact. The first two fields of every non-comment input row are read as
vertex IDs; later fields (for example, weights) are ignored.

The conversion makes a simple undirected graph: self-loops are removed,
directions are ignored, parallel/reciprocal edges are deduplicated, endpoints
are ordered, and adjacency rows are sorted. Vertex IDs must already form the
range selected by GRAPH_VERTEX_BASE and NUM_VERTICES.

Environment:
  GRAPH_VERTEX_BASE  Input ID base, 0 or 1 (default: 0)
  GRAPH_FORCE        Set to 1 to replace existing outputs (default: 0)

NUM_VERTICES preserves isolated vertices. If it is omitted, the value is
inferred as max(normalized vertex ID) + 1; an edgeless input therefore requires
NUM_VERTICES. Lines that are blank or whose first nonblank character is # or %
are ignored.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  usage
  exit 0
fi
if [[ $# -lt 3 || $# -gt 4 ]]; then
  usage >&2
  exit 2
fi

raw_input=$1
adjacency_output=$2
edge_output=$3
declared_vertices=${4:-}
vertex_base=${GRAPH_VERTEX_BASE:-0}
force=${GRAPH_FORCE:-0}

if [[ ! -f $raw_input ]]; then
  echo "Input edge list does not exist: $raw_input" >&2
  exit 2
fi
if [[ $vertex_base != 0 && $vertex_base != 1 ]]; then
  echo "GRAPH_VERTEX_BASE must be 0 or 1." >&2
  exit 2
fi
if [[ $force != 0 && $force != 1 ]]; then
  echo "GRAPH_FORCE must be 0 or 1." >&2
  exit 2
fi
if [[ -n $declared_vertices && ! $declared_vertices =~ ^[0-9]+$ ]]; then
  echo "NUM_VERTICES must be a nonnegative integer." >&2
  exit 2
fi
if [[ $adjacency_output == "$edge_output" ||
      $raw_input == "$adjacency_output" ||
      $raw_input == "$edge_output" ]]; then
  echo "Input, adjacency output, and edge output must be distinct paths." >&2
  exit 2
fi
if [[ $force == 0 ]]; then
  if [[ -e $adjacency_output ]]; then
    echo "Refusing to replace existing output: $adjacency_output" >&2
    exit 2
  fi
  if [[ -e $edge_output ]]; then
    echo "Refusing to replace existing output: $edge_output" >&2
    exit 2
  fi
fi
for tool in awk sort wc mktemp mkdir dirname cp; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required tool is unavailable: $tool" >&2
    exit 2
  fi
done

temporary_root=${TMPDIR:-/tmp}
temporary_dir=$(mktemp -d "$temporary_root/adaptivebk-prepare.XXXXXX")
cleanup() {
  rm -rf -- "$temporary_dir"
}
trap cleanup EXIT

unsorted_edges="$temporary_dir/edges.unsorted"
normalized_edges="$temporary_dir/edges.normalized"
directed_edges="$temporary_dir/edges.directed"
adjacency_file="$temporary_dir/graph.adjacency"

awk -v base="$vertex_base" '
  /^[[:space:]]*$/ || /^[[:space:]]*[#%]/ { next }
  {
    if (NF < 2 || $1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/) {
      printf "invalid edge row at input line %d: %s\n", NR, $0 > "/dev/stderr"
      bad = 1
      next
    }
    u = $1 - base
    v = $2 - base
    if (u < 0 || v < 0) {
      printf "vertex ID below GRAPH_VERTEX_BASE at input line %d\n", NR > "/dev/stderr"
      bad = 1
      next
    }
    if (u == v)
      next
    if (u < v)
      print u, v
    else
      print v, u
  }
  END { if (bad) exit 1 }
' "$raw_input" >"$unsorted_edges"

LC_ALL=C sort -n -k1,1 -k2,2 -u "$unsorted_edges" >"$normalized_edges"

maximum_vertex=$(awk '
  {
    if (!seen || $1 > maximum) maximum = $1
    if ($2 > maximum) maximum = $2
    seen = 1
  }
  END { if (seen) print maximum }
' "$normalized_edges")

if [[ -n $declared_vertices ]]; then
  num_vertices=$declared_vertices
else
  if [[ -z $maximum_vertex ]]; then
    echo "Cannot infer NUM_VERTICES from an edgeless input." >&2
    exit 2
  fi
  num_vertices=$((maximum_vertex + 1))
fi
if [[ -n $maximum_vertex && $maximum_vertex -ge $num_vertices ]]; then
  echo "Normalized vertex ID $maximum_vertex is outside 0..$((num_vertices - 1))." >&2
  exit 2
fi

num_edges=$(awk 'END { print NR + 0 }' "$normalized_edges")
awk '{ print $1, $2; print $2, $1 }' "$normalized_edges" |
  LC_ALL=C sort -n -k1,1 -k2,2 >"$directed_edges"

{
  printf '%s %s\n' "$num_vertices" "$num_edges"
  awk -v n="$num_vertices" '
    BEGIN {
      current = 0
      if (n > 0) line = "0"
    }
    {
      while (current < $1) {
        print line
        current++
        line = sprintf("%d", current)
      }
      if ($1 != current) {
        print "internal error while constructing adjacency rows" > "/dev/stderr"
        failed = 1
        exit 1
      }
      line = line " " $2
    }
    END {
      if (failed) exit 1
      while (current < n) {
        print line
        current++
        if (current < n) line = sprintf("%d", current)
      }
    }
  ' "$directed_edges"
} >"$adjacency_file"

mkdir -p "$(dirname -- "$adjacency_output")" "$(dirname -- "$edge_output")"
cp "$adjacency_file" "$adjacency_output"
cp "$normalized_edges" "$edge_output"

echo "PREPARED vertices=$num_vertices edges=$num_edges adjacency=$adjacency_output edge_list=$edge_output"
