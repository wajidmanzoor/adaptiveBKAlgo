#!/bin/bash

# Usage: ./run_benchmark.sh <data_dir> [logs_dir]
# Defaults: data_dir=/data/labdata/wajid/AdjacencyList, logs_dir=./logs

DATA_DIR="${1:-/data/labdata/wajid/AdjacencyList}"
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

mkdir -p "$LOGS_DIR/pivot"
mkdir -p "$LOGS_DIR/reorder"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTFILE="$LOGS_DIR/benchmark_${TIMESTAMP}.tsv"

printf "Graph\tCliques\tMaxSize\tPivot Checks\tReorder Checks\tPivot Time (ms)\tReorder Time (ms)\n" > "$OUTFILE"

for GRAPH in "$DATA_DIR"/*.txt; do
  [ -f "$GRAPH" ] || continue
  NAME=$(basename "$GRAPH")
  STEM="${NAME%.txt}"

  "$BIN" "$GRAPH" 2 > "$LOGS_DIR/pivot/${STEM}.log" 2>&1
  PIVOT_OUT=$(cat "$LOGS_DIR/pivot/${STEM}.log")
  PIVOT_CLIQUES=$(echo "$PIVOT_OUT" | awk '/Total Maximal Cliques Found:/{print $NF}')
  PIVOT_MAXSIZE=$(echo "$PIVOT_OUT" | awk '/Maximum Clique Size:/{print $NF}')
  PIVOT_CHECKS=$(echo  "$PIVOT_OUT" | awk '/Total Vertex-Set Checks:/{print $NF}')
  PIVOT_TIME=$(echo    "$PIVOT_OUT" | awk '/^Time:/{print $2}')

  "$BIN" "$GRAPH" 5 > "$LOGS_DIR/reorder/${STEM}.log" 2>&1
  REORDER_LINE=$(grep "^ReorderNew:" "$LOGS_DIR/reorder/${STEM}.log")
  REORDER_CLIQUES=$(echo "$REORDER_LINE" | awk -F'cliques=' '{print $2}' | awk '{print $1}')
  REORDER_MAXSIZE=$(echo "$REORDER_LINE" | awk -F'maxSize='  '{print $2}' | awk '{print $1}')
  REORDER_CHECKS=$(echo  "$REORDER_LINE" | awk -F'checks='   '{print $2}' | awk '{print $1}')
  REORDER_TIME=$(echo    "$REORDER_LINE" | awk -F'time='     '{print $2}' | awk '{print $1}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$NAME" "$PIVOT_CLIQUES" "$PIVOT_MAXSIZE" \
    "$PIVOT_CHECKS" "$REORDER_CHECKS" \
    "$PIVOT_TIME" "$REORDER_TIME" >> "$OUTFILE"

  echo "  $NAME: cliques=$PIVOT_CLIQUES  pivot_checks=$PIVOT_CHECKS  reorder_checks=$REORDER_CHECKS  pivot=${PIVOT_TIME}ms  reorder=${REORDER_TIME}ms"
done

echo ""
echo "Results saved to: $OUTFILE"
echo "Per-graph logs:   $LOGS_DIR/pivot/  and  $LOGS_DIR/reorder/"
