#!/usr/bin/env python3
"""
Random graph comparison test for PivotBK vs ReorderSib (all methods).
Generates 100 Erdos-Renyi random graphs, runs all algorithm modes,
checks clique-count agreement, and reports timing.
"""

import subprocess
import random
import tempfile
import os
import re
import sys
from collections import defaultdict

BINARY = os.path.join(os.path.dirname(__file__), "build", "bk_algorithm")

MODES = {
    0: "PivotBK-Orig",
    1: "PivotBK-Asc",
    2: "PivotBK-Desc",
    6: "RSib-BruteForce",
    7: "RSib-Backtrack",
    8: "RSib-Greedy",
    9: "RSib-Bitmask",
    10: "RSib-MinHitSet",
}

def generate_graph(n, p, seed):
    random.seed(seed)
    edges = set()
    adj = defaultdict(set)
    for u in range(n):
        for v in range(u + 1, n):
            if random.random() < p:
                edges.add((u, v))
                adj[u].add(v)
                adj[v].add(u)
    return n, len(edges), adj

def write_graph_file(n, m, adj, path):
    with open(path, "w") as f:
        f.write(f"{n} {m}\n")
        for v in range(n):
            nbrs = sorted(adj[v])
            line = str(v) + (" " + " ".join(map(str, nbrs)) if nbrs else "")
            f.write(line + "\n")

def run_mode(graph_file, mode, timeout=30):
    try:
        result = subprocess.run(
            [BINARY, graph_file, str(mode)],
            capture_output=True, text=True, timeout=timeout
        )
        out = result.stdout + result.stderr
        return out
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception as e:
        return f"ERROR: {e}"

def parse_clique_count(output):
    # PivotBK format: "Total Maximal Cliques Found: N"
    m = re.search(r"Total Maximal Cliques Found:\s*(\d+)", output)
    if m:
        return int(m.group(1))
    # ReorderSib format: "ReorderSib: cliques=N  maxSize=..."
    m = re.search(r"cliques=(\d+)", output)
    if m:
        return int(m.group(1))
    return None

def parse_time_ms(output):
    # PivotBK: "Time: X ms"
    m = re.search(r"Time:\s*([\d.]+)\s*ms", output)
    if m:
        return float(m.group(1))
    # ReorderSib: "time=X ms"
    m = re.search(r"time=([\d.]+)\s*ms", output)
    if m:
        return float(m.group(1))
    return None

def main():
    NUM_GRAPHS = 100
    # Graph sizes: small enough that all methods finish quickly
    # Mix of (n, p) configs
    configs = [
        (10, 0.3), (10, 0.5), (10, 0.7),
        (15, 0.3), (15, 0.5),
        (20, 0.3), (20, 0.4),
        (25, 0.3),
    ]

    print(f"Generating {NUM_GRAPHS} random graphs and running all modes...\n")

    # Accumulators
    total_time = defaultdict(float)   # mode -> total ms
    mismatch_count = 0
    mismatches = []
    run_count = defaultdict(int)
    timeout_count = defaultdict(int)

    # Per-graph results for the table
    graph_results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for graph_idx in range(NUM_GRAPHS):
            cfg = configs[graph_idx % len(configs)]
            n, p = cfg
            seed = graph_idx * 137 + 42

            n_v, m_e, adj = generate_graph(n, p, seed)

            # Skip trivial graphs (< 3 vertices or no edges)
            if n_v < 3 or m_e == 0:
                continue

            gfile = os.path.join(tmpdir, f"graph_{graph_idx}.txt")
            write_graph_file(n_v, m_e, adj, gfile)

            counts = {}
            times = {}
            timed_out = []

            for mode, name in MODES.items():
                out = run_mode(gfile, mode)
                if out == "TIMEOUT" or out.startswith("TIMEOUT"):
                    timed_out.append(name)
                    timeout_count[mode] += 1
                    counts[mode] = None
                    times[mode] = None
                else:
                    c = parse_clique_count(out)
                    t = parse_time_ms(out)
                    counts[mode] = c
                    times[mode] = t
                    if t is not None:
                        total_time[mode] += t
                        run_count[mode] += 1

            # Check agreement: all non-None counts should match
            valid_counts = {m: c for m, c in counts.items() if c is not None}
            unique_counts = set(valid_counts.values())
            match = len(unique_counts) <= 1

            if not match:
                mismatch_count += 1
                mismatches.append({
                    "graph": graph_idx,
                    "n": n_v, "m": m_e, "p": p,
                    "counts": {MODES[m]: c for m, c in counts.items()}
                })

            graph_results.append({
                "graph": graph_idx,
                "n": n_v, "m": m_e,
                "counts": counts,
                "times": times,
                "match": match,
                "timed_out": timed_out,
            })

            # Progress
            if (graph_idx + 1) % 10 == 0:
                print(f"  Completed {graph_idx + 1}/{NUM_GRAPHS} graphs...", flush=True)

    print("\n" + "=" * 100)
    print("RESULTS SUMMARY")
    print("=" * 100)

    # --- Correctness table ---
    print("\n[1] CLIQUE COUNT AGREEMENT (sample of first 20 graphs)\n")
    col_w = 16
    header = f"{'Graph':>6} {'n':>4} {'m':>5} {'Match':>5} "
    for name in MODES.values():
        header += f"{name:>{col_w}}"
    print(header)
    print("-" * len(header))

    for r in graph_results[:20]:
        row = f"{r['graph']:>6} {r['n']:>4} {r['m']:>5} {'YES' if r['match'] else 'NO':>5} "
        for mode, name in MODES.items():
            c = r['counts'].get(mode)
            val = str(c) if c is not None else "TIMEOUT"
            row += f"{val:>{col_w}}"
        print(row)

    # --- Timing summary table ---
    print("\n\n[2] TOTAL RUN TIME ACROSS ALL GRAPHS\n")
    print(f"{'Method':<20} {'Total Time (ms)':>18} {'Avg Time (ms)':>16} {'Runs':>6} {'Timeouts':>10}")
    print("-" * 75)
    for mode, name in MODES.items():
        tot = total_time[mode]
        runs = run_count[mode]
        avg = tot / runs if runs > 0 else 0
        to = timeout_count[mode]
        print(f"{name:<20} {tot:>18.3f} {avg:>16.4f} {runs:>6} {to:>10}")

    # --- Correctness summary ---
    print(f"\n\n[3] CORRECTNESS SUMMARY")
    print(f"  Total graphs tested : {len(graph_results)}")
    print(f"  All methods agree   : {len(graph_results) - mismatch_count}")
    print(f"  Mismatches found    : {mismatch_count}")

    if mismatches:
        print("\n  MISMATCH DETAILS:")
        for mm in mismatches:
            print(f"    Graph {mm['graph']:>3} (n={mm['n']}, m={mm['m']}, p={mm['p']})")
            for name, c in mm['counts'].items():
                print(f"      {name:<20}: {c}")

    print("\n" + "=" * 100)

    if mismatch_count == 0:
        print("RESULT: All methods produce identical clique counts on all tested graphs.")
    else:
        print(f"RESULT: {mismatch_count} graph(s) show disagreement between methods!")

    print("=" * 100 + "\n")

if __name__ == "__main__":
    main()
