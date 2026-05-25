#!/usr/bin/env bash

if [[ -z "${SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BK_SRC_DIR="$REPO_ROOT/orginal"
BK_BUILD_DIR="$BK_SRC_DIR/build"
BK_BIN="$BK_BUILD_DIR/bk_algorithm"
DEFAULT_DATA_DIR="$REPO_ROOT/data"

HBBMC_ROOT="$REPO_ROOT/compare/HBBMC"
HBBMC_SRC_DIR="$HBBMC_ROOT/src"
HBBMC_BUILD_DIR="$HBBMC_SRC_DIR/build"
HBBMC_BIN="$HBBMC_BUILD_DIR/MCE"

job_count() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.logicalcpu 2>/dev/null || echo 4
    else
        echo 4
    fi
}

setup_timeout_runner() {
    local seconds="${1:-120}"
    if command -v gtimeout >/dev/null 2>&1; then
        RUN() { gtimeout "$seconds" "$@"; }
    elif command -v timeout >/dev/null 2>&1; then
        RUN() { timeout "$seconds" "$@"; }
    else
        RUN() { "$@"; }
    fi
}

build_bk_algorithm() {
    cmake -S "$BK_SRC_DIR" -B "$BK_BUILD_DIR" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 \
      && cmake --build "$BK_BUILD_DIR" -j"$(job_count)" >/dev/null 2>&1
}

build_hbbmc() {
    cmake -S "$HBBMC_SRC_DIR" -B "$HBBMC_BUILD_DIR" -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 \
      && cmake --build "$HBBMC_BUILD_DIR" -j"$(job_count)" >/dev/null 2>&1
}
