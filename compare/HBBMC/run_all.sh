#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MCE="$SCRIPT_DIR/src/build/MCE"

usage() {
    echo "Usage: $0 <dataset_dir> <t>"
    echo "  dataset_dir  directory containing graph dataset files"
    echo "  t            early-termination plex value (e.g. 2 or 3)"
    exit 1
}

[[ $# -lt 2 ]] && usage

DATASET_DIR="$1"
T="$2"

if [[ ! -d "$DATASET_DIR" ]]; then
    echo "Error: '$DATASET_DIR' is not a directory." >&2
    exit 1
fi

if [[ ! -f "$MCE" ]]; then
    echo "Error: MCE executable not found at '$MCE'." >&2
    echo "Build it first:" >&2
    echo "  cd src && mkdir -p build && cd build && cmake .. && make" >&2
    exit 1
fi

if [[ ! -x "$MCE" ]]; then
    chmod +x "$MCE"
fi

# Progress bar writes directly to /dev/tty so it shows in terminal
# even when stdout/stderr are redirected to a file via nohup
if [ -w /dev/tty ]; then
    PROGRESS_OUT="/dev/tty"
else
    PROGRESS_OUT="/dev/null"
fi

progress_bar() {
    local current="$1"
    local total="$2"
    local label="$3"
    local width=40

    local percent=$(( current * 100 / total ))
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))

    local bar=""
    for ((i=0; i<filled; i++)); do bar+="#"; done
    for ((i=0; i<empty; i++)); do bar+="-"; done

    printf "\r[%s] %3d%%  %d/%d  %s" "$bar" "$percent" "$current" "$total" "$label" > "$PROGRESS_OUT"

    if [ "$current" -eq "$total" ]; then
        printf "\n" > "$PROGRESS_OUT"
    fi
}

mapfile -t datasets < <(find "$DATASET_DIR" -maxdepth 1 -type f | sort)

N=${#datasets[@]}
if [[ "$N" -eq 0 ]]; then
    echo "No dataset files found in '$DATASET_DIR'." >&2
    exit 1
fi

echo "MCE executable : $MCE"
echo "Dataset dir    : $DATASET_DIR"
echo "t value        : $T"
echo "Datasets found : $N"
echo "----------------------------------------"

completed=0

for dataset in "${datasets[@]}"; do
    name="$(basename "$dataset")"

    progress_bar "$completed" "$N" "Running: $name"

    echo ""
    echo "==> [$((completed + 1))/$N] $name"
    "$MCE" "$dataset" "$T"

    completed=$(( completed + 1 ))
    progress_bar "$completed" "$N" "Done: $name"
done

echo ""
echo "----------------------------------------"
echo "All datasets complete."
