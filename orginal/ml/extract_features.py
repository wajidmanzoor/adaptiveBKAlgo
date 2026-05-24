#!/usr/bin/env python3
"""
extract_features.py — extract structural graph features for ML-based rule selection.

Usage:
    python extract_features.py <graph_dir> [output_csv]

Graph file format (matches C++ loader):
    Line 1:  n m
    Lines 2+: v u1 u2 ...  (vertex v followed by its neighbors)

Output:
    CSV saved to ml/data/graph_features.csv (or the path given as second arg).

Features extracted (no clique or triangle information):
    Basic:       n, m, density, avg/max/min/std degree, skewness, gini, concentration
    Core:        degeneracy, avg_core, std_core, num_core_levels,
                 frac_in_max_core, degeneracy/avg_degree
    Forward deg: max_forward_degree, std_forward_degree  (under degeneracy order)
"""

import sys
import csv
import math
import numpy as np
from pathlib import Path


# ── Graph parsing ─────────────────────────────────────────────────────────────

def parse_graph(filepath):
    """
    Parse the custom adjacency-list format.
    Returns (n, adj) where adj[v] is a set of neighbour vertex IDs.
    """
    with open(filepath) as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        raise ValueError("Empty file")

    first = lines[0].split()
    n, m = int(first[0]), int(first[1])

    adj = [set() for _ in range(n)]
    for line in lines[1:]:
        tokens = list(map(int, line.split()))
        if len(tokens) < 1:
            continue
        v = tokens[0]
        for u in tokens[1:]:
            if u != v and 0 <= u < n and 0 <= v < n:
                adj[v].add(u)
                adj[u].add(v)

    return n, adj


# ── Core decomposition (Batagelj-Zaversnik O(m)) ─────────────────────────────

def core_decomposition(n, adj):
    """
    Compute core number for every vertex using the O(m) bin-sort peeling
    algorithm.  Returns a list core[v] = core number of vertex v.
    """
    deg = [len(adj[v]) for v in range(n)]
    if n == 0:
        return []

    max_deg = max(deg)

    # Build bin start positions (sorted by degree)
    bin_count = [0] * (max_deg + 1)
    for d in deg:
        bin_count[d] += 1

    bin_start = [0] * (max_deg + 1)
    acc = 0
    for d in range(max_deg + 1):
        bin_start[d] = acc
        acc += bin_count[d]

    # Place vertices into sorted order
    starts = list(bin_start)
    pos   = [0] * n   # pos[v]  = index of v in `order`
    order = [0] * n   # order[i] = vertex at position i
    for v in range(n):
        pos[v]         = starts[deg[v]]
        order[pos[v]]  = v
        starts[deg[v]] += 1

    core = list(deg)

    for i in range(n):
        v = order[i]
        for u in adj[v]:
            if core[u] > core[v]:
                du = core[u]
                pu = pos[u]
                pw = bin_start[du]   # first position of bin du
                w  = order[pw]
                if u != w:           # swap u and w inside the bin
                    pos[u], pos[w]     = pw, pu
                    order[pu], order[pw] = w, u
                bin_start[du] += 1   # shrink bin du from the left
                core[u]       -= 1   # u's effective degree drops

    return core


# ── Statistic helpers ─────────────────────────────────────────────────────────

def gini_coefficient(values):
    """Gini coefficient of a non-negative sequence (0 = perfectly equal)."""
    arr = sorted(float(v) for v in values)
    n   = len(arr)
    if n == 0 or sum(arr) == 0:
        return 0.0
    total = sum(arr)
    cum   = sum((2 * (i + 1) - n - 1) * v for i, v in enumerate(arr))
    return cum / (n * total)


def skewness(values):
    """Fisher-Pearson standardised skewness."""
    arr = np.asarray(values, dtype=float)
    std = arr.std()
    if std == 0:
        return 0.0
    return float(((arr - arr.mean()) ** 3).mean() / std ** 3)


# ── Feature extraction ────────────────────────────────────────────────────────

def extract_features(filepath):
    n, adj = parse_graph(filepath)

    if n == 0:
        return None

    degrees = [len(adj[v]) for v in range(n)]
    m       = sum(degrees) // 2

    # ── Basic degree statistics ───────────────────────────────────────────────
    avg_deg     = sum(degrees) / n
    max_deg     = max(degrees)
    min_deg     = min(degrees)
    std_deg     = float(np.std(degrees))
    deg_skew    = skewness(degrees)
    deg_gini    = gini_coefficient(degrees)
    concentration = max_deg / avg_deg if avg_deg > 0 else 0.0
    density     = (2 * m) / (n * (n - 1)) if n > 1 else 0.0

    # ── Core decomposition ────────────────────────────────────────────────────
    core            = core_decomposition(n, adj)
    degeneracy      = max(core)
    avg_core        = float(np.mean(core))
    std_core        = float(np.std(core))
    num_core_levels = len(set(core))
    frac_max_core   = sum(1 for c in core if c == degeneracy) / n
    deg_over_avg    = degeneracy / avg_deg if avg_deg > 0 else 0.0

    # ── Forward degrees under degeneracy order ────────────────────────────────
    # Degeneracy order: sort vertices by core number ascending (peel order).
    # N+(v) = neighbours of v that appear later in this order.
    order_by_core = sorted(range(n), key=lambda v: core[v])
    position      = [0] * n
    for i, v in enumerate(order_by_core):
        position[v] = i

    forward_degs = [
        sum(1 for u in adj[v] if position[u] > position[v])
        for v in range(n)
    ]
    max_fwd_deg = max(forward_degs)
    std_fwd_deg = float(np.std(forward_degs))

    return {
        "graph":                   Path(filepath).stem,
        "n":                       n,
        "m":                       m,
        "density":                 round(density, 8),
        "avg_degree":              round(avg_deg, 4),
        "max_degree":              max_deg,
        "min_degree":              min_deg,
        "std_degree":              round(std_deg, 4),
        "degree_skewness":         round(deg_skew, 4),
        "degree_gini":             round(deg_gini, 4),
        "max_over_avg_degree":     round(concentration, 4),
        "degeneracy":              degeneracy,
        "avg_core":                round(avg_core, 4),
        "std_core":                round(std_core, 4),
        "num_core_levels":         num_core_levels,
        "frac_in_max_core":        round(frac_max_core, 4),
        "degeneracy_over_avg_deg": round(deg_over_avg, 4),
        "max_forward_degree":      max_fwd_deg,
        "std_forward_degree":      round(std_fwd_deg, 4),
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_features.py <graph_dir> [output_csv]")
        sys.exit(1)

    graph_dir  = Path(sys.argv[1])
    script_dir = Path(__file__).parent
    out_dir    = script_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = Path(sys.argv[2]) if len(sys.argv) >= 3 else out_dir / "graph_features.csv"

    graph_files = sorted(graph_dir.glob("*.txt"))
    if not graph_files:
        print(f"No .txt files found in {graph_dir}")
        sys.exit(1)

    print(f"Found {len(graph_files)} graphs in {graph_dir}")

    rows = []
    for gf in graph_files:
        print(f"  {gf.name:<40}", end=" ", flush=True)
        try:
            feat = extract_features(gf)
            if feat is None:
                print("SKIP (empty graph)")
            else:
                rows.append(feat)
                print(f"n={feat['n']:>6}  m={feat['m']:>7}  δ={feat['degeneracy']}")
        except Exception as exc:
            print(f"ERROR — {exc}")

    if not rows:
        print("No features extracted.")
        sys.exit(1)

    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nExtracted {len(rows)} graphs  →  {out_path}")


if __name__ == "__main__":
    main()
