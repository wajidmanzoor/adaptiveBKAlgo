#!/usr/bin/env bash
# compare_with_ours.sh
#
# Compares HBBMC (MCE) vs PivotBK ASC on all graphs in a directory.
# Usage: ./compare_with_ours.sh <data_dir> <our_bin>
#
# Expects:
#   <data_dir>  directory of edge-list (.clean) files for HBBMC
#   <our_bin>   path to our bk_algorithm binary
#
# Both binaries must already be compiled before running this script.
#
# To convert our adjacency-list files to edge-list format first, run:
#   python3 convert_to_edgelist.py <adj_list_dir>
# which outputs to ./edgeListCleaned/

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MCE="$SCRIPT_DIR/src/build/MCE"

usage() {
    echo "Usage: $0 <edgelist_dir> <our_bk_binary>"
    echo ""
    echo "  edgelist_dir   directory of .clean edge-list files (input to HBBMC)"
    echo "  our_bk_binary  path to our compiled bk_algorithm binary"
    exit 1
}

[[ $# -lt 2 ]] && usage

EDIR="$(cd "$1" && pwd)"
OUR_BIN="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"

if [[ ! -d "$EDIR" ]]; then
    echo "Error: '$EDIR' is not a directory." >&2; exit 1
fi
if [[ ! -f "$MCE" ]]; then
    echo "Error: MCE not found at '$MCE'. Build it first:" >&2
    echo "  cd src && mkdir -p build && cd build && cmake .. && make" >&2
    exit 1
fi
if [[ ! -f "$OUR_BIN" ]]; then
    echo "Error: our binary not found at '$OUR_BIN'." >&2; exit 1
fi

if   command -v gtimeout &>/dev/null; then RUN() { gtimeout 120 "$@"; }
elif command -v  timeout &>/dev/null; then RUN() { timeout  120 "$@"; }
else                                        RUN() { "$@"; }
fi

get_hbbmc_count() {
    printf '%s\n' "$1" | grep -oE '#mc = [0-9]+' | grep -oE '[0-9]+$' | head -1
}

get_our_count() {
    local out="$1"
    local c
    c=$(printf '%s\n' "$out" | grep -oE 'Cliques Found: [0-9]+' | grep -oE '[0-9]+$') || true
    printf '%s' "${c:-0}"
}

mapfile -t FILES < <(find "$EDIR" -maxdepth 1 -type f | sort)
N=${#FILES[@]}
[[ "$N" -eq 0 ]] && { echo "No files found in '$EDIR'." >&2; exit 1; }

PASS=0; FAIL=0; TIMEOUT=0

printf "%-35s  %10s  %10s  %6s\n" "Graph" "HBBMC" "Ours" "Status"
printf "%-35s  %10s  %10s  %6s\n" "-----------------------------------" "----------" "----------" "------"

for f in "${FILES[@]}"; do
    name="$(basename "$f")"

    hbbmc_out=$(RUN "$MCE" "$f" 1 2>/dev/null) && rc_h=0 || rc_h=$?
    if [[ "$rc_h" -eq 124 ]]; then
        printf "%-35s  %10s  %10s  %6s\n" "$name" "TIMEOUT" "-" "SKIP"
        ((TIMEOUT++)) || true; continue
    fi
    hbbmc_c=$(get_hbbmc_count "$hbbmc_out")

    # Our algorithm needs the adjacency-list format; derive path by stripping .clean
    adj_name="${name%.clean}"
    adj_file="$SCRIPT_DIR/../../data/$adj_name"
    if [[ ! -f "$adj_file" ]]; then
        printf "%-35s  %10s  %10s  %6s\n" "$name" "${hbbmc_c:-?}" "MISSING" "SKIP"
        continue
    fi

    our_out=$(RUN "$OUR_BIN" "$adj_file" 0 1 0 2>/dev/null) && rc_o=0 || rc_o=$?
    if [[ "$rc_o" -eq 124 ]]; then
        printf "%-35s  %10s  %10s  %6s\n" "$name" "${hbbmc_c:-?}" "TIMEOUT" "SKIP"
        ((TIMEOUT++)) || true; continue
    fi
    our_c=$(get_our_count "$our_out")

    if [[ "${hbbmc_c:-X}" == "${our_c:-Y}" ]]; then
        status="OK"; ((PASS++)) || true
    else
        status="FAIL"; ((FAIL++)) || true
    fi

    printf "%-35s  %10s  %10s  %6s\n" "$name" "${hbbmc_c:-?}" "${our_c:-?}" "$status"
done

echo ""
printf "Results: %d OK  |  %d FAIL  |  %d TIMEOUT  (out of %d)\n" \
    "$PASS" "$FAIL" "$TIMEOUT" "$N"
