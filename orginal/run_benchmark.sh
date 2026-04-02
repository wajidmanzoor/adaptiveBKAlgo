#!/bin/bash

# Usage: ./run_benchmark.sh <data_dir> [logs_dir]
# Defaults: data_dir=./data, logs_dir=./logs

DATA_DIR="${1:- /data/labdata/wajid/AdjacencyList}"
LOGS_DIR="${2:-./logs}"
BIN="./bk_algorithm"

if [ ! -f "$BIN" ]; then
  echo "Error: $BIN not found. Run 'make' first."
  exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
  echo "Error: data directory '$DATA_DIR' not found."
  exit 1
fi

mkdir -p "$LOGS_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTFILE="$LOGS_DIR/benchmark_${TIMESTAMP}.tsv"

printf "Graph\tCliques\tMaxSize\tPivot Checks\tReorder Checks\tPivot Time (ms)\tReorder Time (ms)\n" > "$OUTFILE"

for GRAPH in "$DATA_DIR"/*.txt; do
  [ -f "$GRAPH" ] || continue
  NAME=$(basename "$GRAPH")

  PIVOT_OUT=$("$BIN" "$GRAPH" 2 2>&1)
  PIVOT_CLIQUES=$(echo "$PIVOT_OUT" | grep -oP 'Total Maximal Cliques Found: \K[0-9]+')
  PIVOT_MAXSIZE=$(echo "$PIVOT_OUT" | grep -oP 'Maximum Clique Size: \K[0-9]+')
  PIVOT_CHECKS=$(echo  "$PIVOT_OUT" | grep -oP 'Total Vertex-Set Checks: \K[0-9]+')
  PIVOT_TIME=$(echo    "$PIVOT_OUT" | grep -oP 'Time: \K[0-9]+\.?[0-9]*')

  REORDER_OUT=$("$BIN" "$GRAPH" 5 2>&1)
  REORDER_LINE=$(echo "$REORDER_OUT" | grep "^ReorderNew:")
  REORDER_CLIQUES=$(echo "$REORDER_LINE" | grep -oP 'cliques=\K[0-9]+')
  REORDER_MAXSIZE=$(echo "$REORDER_LINE" | grep -oP 'maxSize=\K[0-9]+')
  REORDER_CHECKS=$(echo  "$REORDER_LINE" | grep -oP 'checks=\K[0-9]+')
  REORDER_TIME=$(echo    "$REORDER_LINE" | grep -oP 'time=\K[0-9]+\.?[0-9]*')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$NAME" "$PIVOT_CLIQUES" "$PIVOT_MAXSIZE" \
    "$PIVOT_CHECKS" "$REORDER_CHECKS" \
    "$PIVOT_TIME" "$REORDER_TIME" >> "$OUTFILE"

  echo "  $NAME: cliques=$PIVOT_CLIQUES  pivot_checks=$PIVOT_CHECKS  reorder_checks=$REORDER_CHECKS  pivot=${PIVOT_TIME}ms  reorder=${REORDER_TIME}ms"
done

echo ""
echo "Results saved to: $OUTFILE"
