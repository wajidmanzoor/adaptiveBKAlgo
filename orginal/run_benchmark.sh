#!/bin/bash

# Usage: ./run_benchmark.sh <data_dir> [logs_dir]
# Defaults: data_dir=/data/labdata/wajid/AdjacencyList, logs_dir=./logs
#
# Modes:
#   2 = Pivot BK        (original order)
#   3 = Pivot BK        (ascending degeneracy)
#   8 = Pivot BK        (descending degeneracy)
#   5 = ReorderNew      (original order)
#   6 = ReorderNew      (ascending degeneracy)
#   7 = ReorderNew      (descending degeneracy)

DATA_DIR="${1:-/data/labdata/wajid/AdjacencyList}"
LOGS_DIR="${2:-./logs}"

if [ -f "./build/bk_algorithm" ]; then
  BIN="./build/bk_algorithm"
elif [ -f "./bk_algorithm" ]; then
  BIN="./bk_algorithm"
else
  echo "Error: bk_algorithm not found. Run 'cmake --build build' or 'make' first."
  exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
  echo "Error: data directory '$DATA_DIR' not found."
  exit 1
fi

mkdir -p "$LOGS_DIR/pivot_orig"
mkdir -p "$LOGS_DIR/pivot_asc"
mkdir -p "$LOGS_DIR/pivot_desc"
mkdir -p "$LOGS_DIR/reorder_orig"
mkdir -p "$LOGS_DIR/reorder_asc"
mkdir -p "$LOGS_DIR/reorder_desc"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTFILE="$LOGS_DIR/benchmark_${TIMESTAMP}.tsv"

printf "Graph\tMaxSize\tP.Orig Cliques\tP.Asc Cliques\tP.Desc Cliques\tR.Orig Cliques\tR.Asc Cliques\tR.Desc Cliques\tP.Orig Checks\tP.Asc Checks\tP.Desc Checks\tR.Orig Checks\tR.Asc Checks\tR.Desc Checks\tP.Orig Time\tP.Asc Time\tP.Desc Time\tR.Orig Time\tR.Asc Time\tR.Desc Time\n" > "$OUTFILE"

for GRAPH in "$DATA_DIR"/*.txt; do
  [ -f "$GRAPH" ] || continue
  NAME=$(basename "$GRAPH")
  STEM="${NAME%.txt}"

  "$BIN" "$GRAPH" 2 > "$LOGS_DIR/pivot_orig/${STEM}.log"   2>&1
  "$BIN" "$GRAPH" 3 > "$LOGS_DIR/pivot_asc/${STEM}.log"    2>&1
  "$BIN" "$GRAPH" 8 > "$LOGS_DIR/pivot_desc/${STEM}.log"   2>&1
  "$BIN" "$GRAPH" 5 > "$LOGS_DIR/reorder_orig/${STEM}.log" 2>&1
  "$BIN" "$GRAPH" 6 > "$LOGS_DIR/reorder_asc/${STEM}.log"  2>&1
  "$BIN" "$GRAPH" 7 > "$LOGS_DIR/reorder_desc/${STEM}.log" 2>&1

  pparse() {
    awk '
      /Total Maximal Cliques Found:/ { c=$NF }
      /Maximum Clique Size:/         { s=$NF }
      /Total Vertex-Set Checks:/     { k=$NF }
      /^Time:/                       { t=$2  }
      END { print c, s, k, t }
    ' "$1"
  }
  rparse() {
    grep "^ReorderNew:" "$1" | awk -F'[ =]+' '{
      for(i=1;i<=NF;i++){
        if($i=="cliques")  c=$(i+1)
        if($i=="maxSize")  s=$(i+1)
        if($i=="checks")   k=$(i+1)
        if($i=="time")     t=$(i+1)
      }
      print c, s, k, t
    }'
  }

  P2C=$(pparse "$LOGS_DIR/pivot_orig/${STEM}.log"  | awk '{print $1}')
  PS=$(pparse  "$LOGS_DIR/pivot_orig/${STEM}.log"  | awk '{print $2}')
  P2K=$(pparse "$LOGS_DIR/pivot_orig/${STEM}.log"  | awk '{print $3}')
  P2T=$(pparse "$LOGS_DIR/pivot_orig/${STEM}.log"  | awk '{print $4}')
  P3C=$(pparse "$LOGS_DIR/pivot_asc/${STEM}.log"   | awk '{print $1}')
  P3K=$(pparse "$LOGS_DIR/pivot_asc/${STEM}.log"   | awk '{print $3}')
  P3T=$(pparse "$LOGS_DIR/pivot_asc/${STEM}.log"   | awk '{print $4}')
  P8C=$(pparse "$LOGS_DIR/pivot_desc/${STEM}.log"  | awk '{print $1}')
  P8K=$(pparse "$LOGS_DIR/pivot_desc/${STEM}.log"  | awk '{print $3}')
  P8T=$(pparse "$LOGS_DIR/pivot_desc/${STEM}.log"  | awk '{print $4}')
  R5C=$(rparse "$LOGS_DIR/reorder_orig/${STEM}.log" | awk '{print $1}')
  R5K=$(rparse "$LOGS_DIR/reorder_orig/${STEM}.log" | awk '{print $3}')
  R5T=$(rparse "$LOGS_DIR/reorder_orig/${STEM}.log" | awk '{print $4}')
  R6C=$(rparse "$LOGS_DIR/reorder_asc/${STEM}.log"  | awk '{print $1}')
  R6K=$(rparse "$LOGS_DIR/reorder_asc/${STEM}.log"  | awk '{print $3}')
  R6T=$(rparse "$LOGS_DIR/reorder_asc/${STEM}.log"  | awk '{print $4}')
  R7C=$(rparse "$LOGS_DIR/reorder_desc/${STEM}.log" | awk '{print $1}')
  R7K=$(rparse "$LOGS_DIR/reorder_desc/${STEM}.log" | awk '{print $3}')
  R7T=$(rparse "$LOGS_DIR/reorder_desc/${STEM}.log" | awk '{print $4}')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$NAME" "$PS" \
    "$P2C" "$P3C" "$P8C" "$R5C" "$R6C" "$R7C" \
    "$P2K" "$P3K" "$P8K" "$R5K" "$R6K" "$R7K" \
    "$P2T" "$P3T" "$P8T" "$R5T" "$R6T" "$R7T" >> "$OUTFILE"

  printf "  %-28s maxSz=%-3s  cliques: po=%-6s pa=%-6s pd=%-6s ro=%-6s ra=%-6s rd=%-6s  checks: po=%-6s pa=%-6s pd=%-6s ro=%-6s ra=%-6s rd=%-6s  time(ms): po=%-10s pa=%-10s pd=%-10s ro=%-10s ra=%-10s rd=%s\n" \
    "$NAME" "$PS" \
    "$P2C" "$P3C" "$P8C" "$R5C" "$R6C" "$R7C" \
    "$P2K" "$P3K" "$P8K" "$R5K" "$R6K" "$R7K" \
    "$P2T" "$P3T" "$P8T" "$R5T" "$R6T" "$R7T"
done

echo ""
echo "Results saved to: $OUTFILE"
