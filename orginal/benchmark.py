#!/usr/bin/env python3
"""
Bottleneck benchmark for ReorderSib methods vs PivotBK.
Runs graphs at increasing sizes/densities, repeated N times each,
to expose which method/graph-size combo is the bottleneck.
"""

import subprocess, random, tempfile, os, re, sys
from collections import defaultdict

BINARY = os.path.join(os.path.dirname(__file__), "build", "release_build", "bk_algorithm")

METHODS = {
    0:  "PivotBK",
    6:  "RSib-BruteForce",
    7:  "RSib-Backtrack",
    8:  "RSib-Greedy",
    9:  "RSib-Bitmask",
    10: "RSib-MinHitSet",
}

def gen_graph(n, p, seed):
    random.seed(seed)
    adj = defaultdict(set)
    edges = 0
    for u in range(n):
        for v in range(u+1, n):
            if random.random() < p:
                adj[u].add(v); adj[v].add(u)
                edges += 1
    return n, edges, adj

def write_graph(n, m, adj, path):
    with open(path, "w") as f:
        f.write(f"{n} {m}\n")
        for v in range(n):
            nbrs = sorted(adj[v])
            f.write(str(v) + (" " + " ".join(map(str, nbrs)) if nbrs else "") + "\n")

def run(binary, gfile, mode, timeout=60):
    try:
        r = subprocess.run([binary, gfile, str(mode)],
                           capture_output=True, text=True, timeout=timeout)
        return r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return "TIMEOUT"

def parse(out):
    c = re.search(r"Total Maximal Cliques Found:\s*(\d+)|cliques=(\d+)", out)
    t = re.search(r"Time:\s*([\d.]+)\s*ms|time=([\d.]+)\s*ms", out)
    cliques = int(c.group(1) or c.group(2)) if c else None
    ms = float(t.group(1) or t.group(2)) if t else None
    return cliques, ms

# ── Graph configs: (n, p, label)
CONFIGS = [
    (15,  0.4, "n=15  p=0.4"),
    (20,  0.4, "n=20  p=0.4"),
    (25,  0.4, "n=25  p=0.4"),
    (30,  0.35,"n=30  p=0.35"),
    (35,  0.3, "n=35  p=0.3"),
    (40,  0.3, "n=40  p=0.3"),
    (50,  0.25,"n=50  p=0.25"),
    (60,  0.2, "n=60  p=0.2"),
    (75,  0.15,"n=75  p=0.15"),
    (100, 0.12,"n=100 p=0.12"),
]
REPS = 5   # graphs per config (different seeds)

def main():
    # total_time[config_label][mode] -> list of ms
    times  = defaultdict(lambda: defaultdict(list))
    counts = defaultdict(lambda: defaultdict(list))

    with tempfile.TemporaryDirectory() as tmp:
        for cfg in CONFIGS:
            n, p, label = cfg
            print(f"  {label} ...", flush=True)
            for rep in range(REPS):
                seed = rep * 997 + n * 13
                nv, me, adj = gen_graph(n, p, seed)
                if me == 0:
                    continue
                gfile = os.path.join(tmp, f"g_{n}_{rep}.txt")
                write_graph(nv, me, adj, gfile)

                for mode in METHODS:
                    out = run(BINARY, gfile, mode)
                    if out == "TIMEOUT":
                        times[label][mode].append(None)
                    else:
                        c, t = parse(out)
                        times[label][mode].append(t)
                        counts[label][mode].append(c)

    # ── Table 1: Average time per config per method
    col = 14
    print("\n" + "=" * 100)
    print("AVERAGE TIME (ms) PER GRAPH SIZE — Release build")
    print("=" * 100)
    header = f"{'Graph config':<18}" + "".join(f"{METHODS[m]:>{col}}" for m in METHODS)
    print(header)
    print("-" * len(header))

    for cfg in CONFIGS:
        label = cfg[2]
        row = f"{label:<18}"
        for mode in METHODS:
            ts = [x for x in times[label][mode] if x is not None]
            if not ts:
                row += f"{'TIMEOUT':>{col}}"
            else:
                row += f"{sum(ts)/len(ts):>{col}.3f}"
        print(row)

    # ── Table 2: Slowdown vs PivotBK
    print("\n\nSLOWDOWN RELATIVE TO PivotBK (×)")
    print("-" * len(header))
    header2 = f"{'Graph config':<18}" + "".join(f"{METHODS[m]:>{col}}" for m in METHODS if m != 0)
    print(f"{'Graph config':<18}" + "".join(f"{METHODS[m]:>{col}}" for m in METHODS))
    print("-" * len(header))

    for cfg in CONFIGS:
        label = cfg[2]
        pivot_ts = [x for x in times[label][0] if x is not None]
        pivot_avg = sum(pivot_ts)/len(pivot_ts) if pivot_ts else None
        row = f"{label:<18}"
        for mode in METHODS:
            ts = [x for x in times[label][mode] if x is not None]
            if not ts or pivot_avg is None:
                row += f"{'N/A':>{col}}"
            elif mode == 0:
                row += f"{'1.0×':>{col}}"
            else:
                ratio = (sum(ts)/len(ts)) / pivot_avg
                row += f"{ratio:>{col}.1f}×"
        print(row)

    # ── Table 3: Correctness check (all counts match PivotBK?)
    print("\n\nCORRECTNESS: does clique count match PivotBK? (across all reps)")
    print("-" * len(header))
    for cfg in CONFIGS:
        label = cfg[2]
        row = f"{label:<18}"
        pivot_counts = counts[label][0]
        for mode in METHODS:
            if mode == 0:
                row += f"{'baseline':>{col}}"
                continue
            mc = counts[label][mode]
            if not mc:
                row += f"{'N/A':>{col}}"
                continue
            # compare element-wise where both are available
            mismatches = sum(1 for a, b in zip(pivot_counts, mc)
                             if a is not None and b is not None and a != b)
            total = sum(1 for a, b in zip(pivot_counts, mc)
                        if a is not None and b is not None)
            if mismatches == 0:
                row += f"{'OK':>{col}}"
            else:
                row += f"{'WRONG('+str(mismatches)+')':>{col}}"
        print(row)

    # ── Table 4: Scaling analysis – how fast does each method grow?
    print("\n\nSCALING: ratio of time at n=100 vs n=15 (shows growth rate)")
    print("-" * 60)
    for mode, name in METHODS.items():
        small = [x for x in times[CONFIGS[0][2]][mode] if x is not None]
        large = [x for x in times[CONFIGS[-1][2]][mode] if x is not None]
        if small and large:
            ratio = (sum(large)/len(large)) / (sum(small)/len(small))
            print(f"  {name:<20}: {ratio:>8.1f}×  ({sum(small)/len(small):.3f} ms → {sum(large)/len(large):.3f} ms)")
        else:
            print(f"  {name:<20}: TIMEOUT on large graphs")

    print("=" * 100)

if __name__ == "__main__":
    print(f"Running benchmark with Release binary: {BINARY}\n")
    main()
