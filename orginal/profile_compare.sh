#!/usr/bin/env bash
# profile_compare.sh — side-by-side profiling: Backtrack+ASC vs Optimized+ASC
# Aggregates per-component time across all graphs and reports % of wall time.
# Usage: ./profile_compare.sh [data_dir]

set -uo pipefail
cd "$(dirname "$0")"

BIN="./bk_algorithm"
DATA="${1:-../data}"

# macOS has no built-in 'timeout'
if command -v gtimeout &>/dev/null; then
    RUN() { gtimeout 120 "$@"; }
elif command -v timeout &>/dev/null; then
    RUN() { timeout 120 "$@"; }
else
    RUN() { "$@"; }
fi

# ── Build ─────────────────────────────────────────────────────────────────────
echo "=== Building ==="
make -C "$(dirname "$0")" -s && echo "OK"
echo ""

# ── Discover graphs ───────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do
    GRAPHS+=("$f")
done < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }
echo "Graphs : $N  ($DATA)"
echo ""

# ── Run one method, dump per-graph profiling lines to a temp file ─────────────
run_method() {
    local label="$1" algo="$2" ord="$3" meth="$4" outfile="$5"
    : > "$outfile"
    local done=0
    for g in "${GRAPHS[@]}"; do
        out=$(RUN "$BIN" "$g" "$algo" "$ord" "$meth" 2>/dev/null) || out=""
        # Append all output — Python will extract the profiling lines
        printf '%s\n' "$out" >> "$outfile"
        (( done++ )) || true
        printf "\r  %s : %d/%d graphs" "$label" "$done" "$N" >&2
    done
    printf "\r  %s : done (%d graphs)          \n" "$label" "$N" >&2
}

TMPBT=$(mktemp)
TMPOPT=$(mktemp)
trap 'rm -f "$TMPBT" "$TMPOPT"' EXIT

echo "Running Backtrack+ASC  (algo=1 ord=1 meth=1)..."
run_method "Backtrack+ASC" 1 1 1 "$TMPBT"

echo "Running Optimized+ASC  (algo=1 ord=1 meth=5)..."
run_method "Optimized+ASC" 1 1 5 "$TMPOPT"

# ── Aggregate and display with Python ────────────────────────────────────────
python3 - "$TMPBT" "$TMPOPT" << 'PYEOF'
import sys, re

# Direct (non-recursive) components — their ms sums are genuine wall-time fractions.
# enumerate/rCall are excluded here because their timers fire on every recursive call,
# so their accumulated ms far exceeds wall time and cannot be used as wall-% fractions.
DIRECT = [
    ("solver",   "solver"),
    ("collect",  "collectCoveringCliques"),
    ("buildhit", "buildHitSets"),
    ("minimal",  "minimalByInclusion"),
    ("common",   "commonExpand"),
    ("dedup",    "seenCliques dedup"),
]

def parse_file(path):
    totals = {k: 0.0 for k, _ in DIRECT}
    totals["wall"] = 0.0
    calls  = {k: 0   for k, _ in DIRECT}

    with open(path) as f:
        for line in f:
            mw = re.search(r'TOTAL.*?([\d.]+) ms', line)
            m  = re.search(r'([\d.]+) ms\s+[\d.]+%\s+calls=(\d+)', line)
            if mw:
                totals["wall"] += float(mw.group(1))
            elif m:
                ms_val = float(m.group(1))
                n_val  = int(m.group(2))
                for key, label in DIRECT:
                    if label.split()[0].lower() in line.lower():
                        totals[key] += ms_val
                        calls[key]  += n_val
                        break
    return totals, calls

bt_t,  bt_c  = parse_file(sys.argv[1])
opt_t, opt_c = parse_file(sys.argv[2])

def pct(v, wall):
    return 100.0 * v / wall if wall > 0 else 0.0

wall_bt  = bt_t["wall"]
wall_opt = opt_t["wall"]

# Sibling-effect subtotal (sum of all direct components)
sib_bt  = sum(bt_t[k]  for k, _ in DIRECT)
sib_opt = sum(opt_t[k] for k, _ in DIRECT)

# Framework = wall - sibling (enumerate + rCall inferred, not double-counted)
fw_bt  = max(0.0, wall_bt  - sib_bt)
fw_opt = max(0.0, wall_opt - sib_opt)

W = 78
print()
print("=" * W)
print(f"  {'Component':<30}  {'Backtrack+ASC':>20}  {'Optimized+ASC':>20}")
print(f"  {'':30}  {'ms':>7} {'%wall':>6} {'calls':>5}  {'ms':>7} {'%wall':>6} {'calls':>5}")
print("-" * W)

# Sibling-effect section
print(f"  {'── Sibling effect ──':<30}")
for key, label in DIRECT:
    bms = bt_t[key];  bc = bt_c[key]
    oms = opt_t[key]; oc = opt_c[key]
    print(f"  {label:<30}  {bms:>7.1f} {pct(bms,wall_bt):>5.1f}% {bc:>5}  "
          f"{oms:>7.1f} {pct(oms,wall_opt):>5.1f}% {oc:>5}")

print(f"  {'  subtotal':<30}  {sib_bt:>7.1f} {pct(sib_bt,wall_bt):>5.1f}%  {'':>5}  "
      f"{sib_opt:>7.1f} {pct(sib_opt,wall_opt):>5.1f}%  {'':>5}")

print()
# Framework section (inferred)
print(f"  {'── Framework (enumerate+rCall) ──':<30}")
print(f"  {'  inferred':<30}  {fw_bt:>7.1f} {pct(fw_bt,wall_bt):>5.1f}%  {'':>5}  "
      f"{fw_opt:>7.1f} {pct(fw_opt,wall_opt):>5.1f}%  {'':>5}")

print("-" * W)
print(f"  {'TOTAL (wall)':<30}  {wall_bt:>7.1f} {'100.0%':>6}  {'':>5}  "
      f"{wall_opt:>7.1f} {'100.0%':>6}  {'':>5}")
print("=" * W)
print()

speedup = wall_bt / wall_opt if wall_opt > 0 else 0
print(f"  Speedup                  : {speedup:.2f}x")
print(f"  Solver share (Backtrack) : {pct(bt_t['solver'], wall_bt):.1f}%  ({bt_t['solver']:.1f} ms)")
print(f"  Solver share (Optimized) : {pct(opt_t['solver'], wall_opt):.1f}%  ({opt_t['solver']:.1f} ms)  →  {bt_t['solver']/opt_t['solver']:.1f}x faster")
print(f"  Framework share (BT)     : {pct(fw_bt, wall_bt):.1f}%  ({fw_bt:.1f} ms)")
print(f"  Framework share (OPT)    : {pct(fw_opt, wall_opt):.1f}%  ({fw_opt:.1f} ms)")
print()
PYEOF
