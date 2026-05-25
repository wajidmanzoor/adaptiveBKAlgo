#!/usr/bin/env bash
# profile_reordersib_big.sh
#
# Build a profiling-enabled ReorderSib binary and run it on the largest graphs
# in a dataset directory, then aggregate the internal cost breakdown.
#
# Usage: ./profile_reordersib_big.sh [data_dir] [top_k]
#   data_dir  graph directory, defaults to repo_root/data
#   top_k     number of largest graphs to profile, defaults to 5

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

DATA_DIR="${1:-$DEFAULT_DATA_DIR}"
TOP_K="${2:-5}"
BIN="$BK_PROFILE_BIN"

setup_timeout_runner 1800

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: '$DATA_DIR' is not a directory." >&2
    exit 1
fi

case "$TOP_K" in
    ''|*[!0-9]*)
        echo "Error: top_k must be a positive integer." >&2
        exit 1
        ;;
esac

build_bk_algorithm_profile && echo "Built profiling bk_algorithm" || {
    echo "Failed to build profiling bk_algorithm" >&2
    exit 1
}
echo ""

TMP_OUT="$(mktemp /private/tmp/reordersib_profile.XXXXXX)"
trap 'rm -f "$TMP_OUT"' EXIT

GRAPHS=()
while IFS= read -r line; do
    [[ -n "$line" ]] && GRAPHS+=("$line")
done < <(python3 - "$DATA_DIR" "$TOP_K" <<'PY'
import sys
from pathlib import Path

data_dir = Path(sys.argv[1])
top_k = int(sys.argv[2])
files = [p for p in data_dir.iterdir() if p.is_file() and not p.name.startswith(".")]
files.sort(key=lambda p: (p.stat().st_size, p.name), reverse=True)
for path in files[:top_k]:
    print(path)
PY
)

N=${#GRAPHS[@]}
if [[ "$N" -eq 0 ]]; then
    echo "No graph files found in $DATA_DIR" >&2
    exit 1
fi

printf "Data dir : %s\n" "$DATA_DIR"
printf "Top K    : %s\n" "$TOP_K"
printf "Graphs   : %d\n\n" "$N"

for g in "${GRAPHS[@]}"; do
    name="$(basename "$g")"
    printf "==> %s\n" "$name"
    out=$(RUN "$BIN" "$g" 1 1 1 4294967295 1 1 1 1 1 1 1 1 2>&1) && rc=0 || rc=$?
    printf '===GRAPH=== %s rc=%s\n' "$g" "$rc" >> "$TMP_OUT"
    printf '%s\n' "$out" >> "$TMP_OUT"
    printf '%s\n\n' "$out"
done

echo ""
python3 - "$TMP_OUT" <<'PY'
import re
import sys
from collections import defaultdict

path = sys.argv[1]

component_re = re.compile(r'^\s{2}(.+?)\s+([0-9.]+) ms\s+([0-9.]+)%\s+calls=(\d+)\s*$')
total_re = re.compile(r'ReorderSib:\s+cliques=\d+\s+dups=\d+\s+maxSize=\d+\s+checks=\d+\s+time=([0-9.]+)\s+ms')
marker_re = re.compile(r'^===GRAPH===\s+(.+)\s+rc=(\d+)$')

graphs = []
current = None

with open(path) as f:
    for raw in f:
        line = raw.rstrip("\n")
        m = marker_re.match(line)
        if m:
            current = {
                "graph": m.group(1),
                "rc": int(m.group(2)),
                "time_ms": None,
                "components": {},
            }
            graphs.append(current)
            continue
        if current is None:
            continue
        m = total_re.search(line)
        if m:
            current["time_ms"] = float(m.group(1))
            continue
        m = component_re.match(line)
        if m:
            current["components"][m.group(1).strip()] = {
                "ms": float(m.group(2)),
                "pct": float(m.group(3)),
                "calls": int(m.group(4)),
            }

ok_graphs = [g for g in graphs if g["rc"] == 0 and g["time_ms"] is not None]
if not ok_graphs:
    print("No successful profiling runs to aggregate.")
    sys.exit(0)

totals = defaultdict(float)
calls = defaultdict(int)
wall = 0.0
for g in ok_graphs:
    wall += g["time_ms"]
    for name, entry in g["components"].items():
        totals[name] += entry["ms"]
        calls[name] += entry["calls"]

print("Aggregate ReorderSib Profile")
print("----------------------------")
print(f"Successful graphs : {len(ok_graphs)} / {len(graphs)}")
print(f"Aggregate wall ms : {wall:.3f}")
print("")
print(f"{'Component':<30} {'Total ms':>12} {'% wall':>10} {'Calls':>10}")
print(f"{'-'*30:<30} {'-'*12:>12} {'-'*10:>10} {'-'*10:>10}")
for name, ms in sorted(totals.items(), key=lambda kv: kv[1], reverse=True):
    pct = 100.0 * ms / wall if wall > 0 else 0.0
    print(f"{name:<30} {ms:>12.3f} {pct:>9.1f}% {calls[name]:>10}")
PY
