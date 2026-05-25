#!/usr/bin/env bash
# check_correctness.sh
#
# Correctness check: compare clique counts across algorithms.
# Reference: PivotBK ASC
# Checked  : ReorderSib Optimized ASC (all rules ON)
#            ReorderSib Optimized ASC (all rules OFF)
#
# Usage: ./check_correctness.sh [data_dir]

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"
cd "$SCRIPT_DIR"

BIN="$BK_BIN"
DATA="${1:-$DEFAULT_DATA_DIR}"

setup_timeout_runner 300

# ── Build ──────────────────────────────────────────────────────────────────────
build_bk_algorithm \
  && echo "Build OK" || { echo "Build FAILED"; exit 1; }

# ── Graph list ─────────────────────────────────────────────────────────────────
GRAPHS=()
while IFS= read -r -d '' f; do GRAPHS+=("$f"); done \
    < <(find "$DATA" -maxdepth 1 -name "*.txt" -print0 | sort -z)
N=${#GRAPHS[@]}
[ "$N" -eq 0 ] && { echo "No graphs in $DATA"; exit 1; }
echo "Graphs : $N   ($DATA)"
echo ""

# ── Helpers ────────────────────────────────────────────────────────────────────
get_cliques() {
    printf '%s\n' "$1" \
        | grep -oE 'Cliques Found: [0-9]+|cliques=[0-9]+' \
        | grep -oE '[0-9]+$' | head -1
}
get_time() {
    printf '%s\n' "$1" \
        | grep -oE 'Time: [0-9]+\.[0-9]+|time=[0-9]+\.[0-9]+' \
        | grep -oE '[0-9.]+$' | head -1
}

# ── Run and compare ────────────────────────────────────────────────────────────
PASS=0
FAIL=0
TIMEOUT=0

FMT="%-35s  %10s  %10s  %10s  %10s  %6s\n"
printf "$FMT" "Graph" "Ref(ms)" "OptON(ms)" "OptOFF(ms)" "Cliques" "Status"
printf "$FMT" "-----------------------------------" "----------" "----------" "----------" "----------" "------"

for g in "${GRAPHS[@]}"; do
    name="$(basename "$g")"

    # Reference: PivotBK ASC
    out_ref=$(RUN "$BIN" "$g" 0 1 0 2>/dev/null) && rc_ref=0 || rc_ref=$?
    ref_c=$(get_cliques "$out_ref")
    ref_t=$(get_time "$out_ref")

    if [ "$rc_ref" -eq 124 ]; then
        printf "%-35s  %10s  %10s  %10s  %10s  %6s\n" \
               "$name" "TIMEOUT" "-" "-" "-" "SKIP"
        ((TIMEOUT++)) || true
        continue
    fi

    # ReorderSib Optimized ASC — all rules ON
    out_on=$(RUN "$BIN" "$g" 1 1 1 4294967295 1 1 1 1 1 1 1 1 2>/dev/null) && rc_on=0 || rc_on=$?
    on_c=$(get_cliques "$out_on")
    on_t=$(get_time "$out_on")

    # ReorderSib Optimized ASC — all rules OFF
    out_off=$(RUN "$BIN" "$g" 1 1 1 4294967295 0 0 0 0 0 0 0 0 2>/dev/null) && rc_off=0 || rc_off=$?
    off_c=$(get_cliques "$out_off")
    off_t=$(get_time "$out_off")

    # Check counts
    ok_on=true;  ok_off=true
    [ "${on_c:-0}"  != "${ref_c:-0}" ] && ok_on=false
    [ "${off_c:-0}" != "${ref_c:-0}" ] && ok_off=false

    if $ok_on && $ok_off; then
        status="OK"
        ((PASS++)) || true
    else
        status="FAIL"
        ((FAIL++)) || true
        [ "$ok_on"  = "false" ] && status="${status}[ON:${on_c}≠${ref_c}]"
        [ "$ok_off" = "false" ] && status="${status}[OFF:${off_c}≠${ref_c}]"
    fi

    printf "%-35s  %10s  %10s  %10s  %10s  %6s\n" \
           "$name" \
           "${ref_t:-?}" \
           "${on_t:-?}" \
           "${off_t:-?}" \
           "${ref_c:-?}" \
           "$status"
done

echo ""
printf "Results: %d OK  |  %d FAIL  |  %d TIMEOUT  (out of %d graphs)\n" \
       "$PASS" "$FAIL" "$TIMEOUT" "$N"
