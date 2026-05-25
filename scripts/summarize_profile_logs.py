#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path


HEADER_RE = re.compile(
    r"ReorderSib:\s+cliques=(\d+)\s+dups=(\d+)\s+maxSize=(\d+)\s+checks=(\d+)\s+time=([0-9.]+)\s+ms"
)
ROW_RE = re.compile(
    r"^\s{2}(.+?)\s+([0-9.]+) ms\s+([0-9.]+)%(\s+calls=(\d+))?\s*$"
)

SKIP_COMPONENTS = {"PROFILED subtotal", "TOTAL (wall)"}


def parse_log(path: Path) -> dict | None:
    text = path.read_text(errors="ignore")
    header = HEADER_RE.search(text)
    if not header:
        return None

    components: dict[str, dict[str, float | int]] = {}
    for line in text.splitlines():
        row = ROW_RE.match(line)
        if not row:
            continue
        name = row.group(1).strip()
        ms = float(row.group(2))
        pct = float(row.group(3))
        calls = int(row.group(5)) if row.group(5) is not None else 0
        components[name] = {"ms": ms, "pct": pct, "calls": calls}

    return {
        "path": path,
        "graph": path.name.removesuffix(".log"),
        "cliques": int(header.group(1)),
        "dups": int(header.group(2)),
        "max_size": int(header.group(3)),
        "checks": int(header.group(4)),
        "wall_ms": float(header.group(5)),
        "components": components,
    }


def summarize(log_dir: Path, top_n: int) -> str:
    logs = []
    for path in sorted(log_dir.glob("*.log")):
        if path.name == "profile_all_data_run.log":
            continue
        parsed = parse_log(path)
        if parsed is not None:
            logs.append(parsed)

    if not logs:
        return f"No per-graph profile logs found in {log_dir}"

    total_wall = sum(item["wall_ms"] for item in logs)
    totals_ms: defaultdict[str, float] = defaultdict(float)
    totals_calls: defaultdict[str, int] = defaultdict(int)
    hottest_count: Counter[str] = Counter()
    hottest_examples: dict[str, tuple[str, float]] = {}

    for item in logs:
        best_name = None
        best_ms = -1.0
        for name, entry in item["components"].items():
            if name in SKIP_COMPONENTS:
                continue
            ms = float(entry["ms"])
            totals_ms[name] += ms
            totals_calls[name] += int(entry["calls"])
            if name != "other / untimed" and ms > best_ms:
                best_name = name
                best_ms = ms
        if best_name is not None:
            hottest_count[best_name] += 1
            prev = hottest_examples.get(best_name)
            if prev is None or best_ms > prev[1]:
                hottest_examples[best_name] = (item["graph"], best_ms)

    lines: list[str] = []
    lines.append("Aggregate ReorderSib Profile Summary")
    lines.append("-----------------------------------")
    lines.append(f"Log dir          : {log_dir}")
    lines.append(f"Parsed graphs    : {len(logs)}")
    lines.append(f"Aggregate wall   : {total_wall:.3f} ms")
    lines.append("")
    lines.append("Top Components By Total Self-Time")
    lines.append(
        f"{'Component':<30} {'Total ms':>12} {'% wall':>9} {'Calls':>12} {'Graphs hot':>11}"
    )
    lines.append(
        f"{'-'*30:<30} {'-'*12:>12} {'-'*9:>9} {'-'*12:>12} {'-'*11:>11}"
    )
    ranked = sorted(totals_ms.items(), key=lambda kv: kv[1], reverse=True)
    for name, ms in ranked[:top_n]:
        pct = 100.0 * ms / total_wall if total_wall else 0.0
        lines.append(
            f"{name:<30} {ms:>12.3f} {pct:>8.1f}% {totals_calls[name]:>12} {hottest_count[name]:>11}"
        )

    lines.append("")
    lines.append("Most Frequent Per-Graph Hotspot")
    lines.append(f"{'Component':<30} {'Graphs':>8} {'Largest single case':>24}")
    lines.append(f"{'-'*30:<30} {'-'*8:>8} {'-'*24:>24}")
    for name, count in hottest_count.most_common(top_n):
        example_graph, example_ms = hottest_examples[name]
        lines.append(
            f"{name:<30} {count:>8} {example_graph + f' ({example_ms:.3f} ms)':>24}"
        )

    lines.append("")
    lines.append("Slowest Graphs")
    lines.append(f"{'Graph':<45} {'Wall ms':>12} {'Cliques':>10} {'Max':>6}")
    lines.append(f"{'-'*45:<45} {'-'*12:>12} {'-'*10:>10} {'-'*6:>6}")
    for item in sorted(logs, key=lambda x: x["wall_ms"], reverse=True)[:top_n]:
        lines.append(
            f"{item['graph']:<45} {item['wall_ms']:>12.3f} {item['cliques']:>10} {item['max_size']:>6}"
        )

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate ReorderSib per-graph profile logs into a corpus-level hotspot summary."
    )
    parser.add_argument(
        "log_dir",
        nargs="?",
        default="scripts/profile_logs",
        help="Directory containing per-graph .log files (default: scripts/profile_logs)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=12,
        help="How many rows to show in each ranked section (default: 12)",
    )
    args = parser.parse_args()

    print(summarize(Path(args.log_dir), args.top))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
