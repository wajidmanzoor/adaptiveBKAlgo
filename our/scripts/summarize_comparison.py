#!/usr/bin/env python3
"""Summarize paired run_external_comparison.sh CSV files."""

from __future__ import annotations

import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def usage() -> None:
    print(
        "usage: summarize_comparison.py OUTPUT_CSV CANDIDATE BASELINE "
        "RUNS_CSV [RUNS_CSV ...]",
        file=sys.stderr,
    )


def main() -> int:
    if len(sys.argv) < 5:
        usage()
        return 2

    output_path = Path(sys.argv[1])
    candidate = sys.argv[2]
    baseline = sys.argv[3]
    input_paths = [Path(argument) for argument in sys.argv[4:]]
    if output_path.exists():
        print(f"refusing to overwrite existing output: {output_path}", file=sys.stderr)
        return 2

    required = {
        "dataset",
        "variant",
        "budget",
        "min_clique_size",
        "status",
        "cliques",
        "runtime_ms",
    }
    all_groups: set[tuple[str, int, int]] = set()
    runtimes: dict[tuple[str, int, int, str], list[float]] = defaultdict(list)
    timeout_counts: dict[tuple[str, int, int, str], int] = defaultdict(int)
    clique_counts: dict[tuple[str, int, int], set[int]] = defaultdict(set)

    for input_path in input_paths:
        with input_path.open(newline="") as source:
            reader = csv.DictReader(source)
            missing = required.difference(reader.fieldnames or ())
            if missing:
                print(
                    f"{input_path}: missing columns: {', '.join(sorted(missing))}",
                    file=sys.stderr,
                )
                return 2
            for row in reader:
                variant = row["variant"]
                if variant not in {candidate, baseline}:
                    continue
                group = (
                    row["dataset"],
                    int(row["budget"]),
                    int(row["min_clique_size"]),
                )
                all_groups.add(group)
                key = (*group, variant)
                if row["status"] == "ok":
                    runtimes[key].append(float(row["runtime_ms"]))
                    clique_counts[group].add(int(row["cliques"]))
                elif row["status"] == "timeout":
                    timeout_counts[key] += 1

    mismatches = {
        group: counts for group, counts in clique_counts.items() if len(counts) > 1
    }
    if mismatches:
        for group, counts in sorted(mismatches.items()):
            print(
                f"clique-count mismatch for {group}: {sorted(counts)}",
                file=sys.stderr,
            )
        return 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    speedups: list[float] = []
    rows: list[list[object]] = []
    for group in sorted(all_groups, key=lambda item: (item[1], item[0], item[2])):
        candidate_key = (*group, candidate)
        baseline_key = (*group, baseline)
        candidate_values = runtimes[candidate_key]
        baseline_values = runtimes[baseline_key]
        candidate_median = (
            statistics.median(candidate_values) if candidate_values else None
        )
        baseline_median = (
            statistics.median(baseline_values) if baseline_values else None
        )
        speedup = (
            baseline_median / candidate_median
            if candidate_median is not None and baseline_median is not None
            else None
        )
        if speedup is not None:
            speedups.append(speedup)
        counts = clique_counts.get(group, set())
        rows.append(
            [
                *group,
                next(iter(counts)) if counts else "",
                len(candidate_values),
                timeout_counts[candidate_key],
                f"{candidate_median:.3f}" if candidate_median is not None else "",
                len(baseline_values),
                timeout_counts[baseline_key],
                f"{baseline_median:.3f}" if baseline_median is not None else "",
                f"{speedup:.4f}" if speedup is not None else "",
            ]
        )

    with output_path.open("x", newline="") as destination:
        writer = csv.writer(destination)
        writer.writerow(
            [
                "dataset",
                "budget",
                "min_clique_size",
                "cliques",
                f"{candidate}_ok_runs",
                f"{candidate}_timeouts",
                f"{candidate}_median_runtime_ms",
                f"{baseline}_ok_runs",
                f"{baseline}_timeouts",
                f"{baseline}_median_runtime_ms",
                f"{candidate}_speedup_vs_{baseline}",
            ]
        )
        writer.writerows(rows)

    wins = sum(speedup > 1.01 for speedup in speedups)
    losses = sum(speedup < 0.99 for speedup in speedups)
    ties = len(speedups) - wins - losses
    geometric_mean = (
        math.exp(statistics.fmean(math.log(speedup) for speedup in speedups))
        if speedups
        else math.nan
    )
    print(
        f"SUMMARY_COMPLETE paired={len(speedups)} wins_gt_1pct={wins} "
        f"losses_gt_1pct={losses} ties_within_1pct={ties} "
        f"geometric_mean_speedup={geometric_mean:.4f} output={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
