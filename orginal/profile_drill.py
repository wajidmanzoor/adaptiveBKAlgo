#!/usr/bin/env python3
"""
Drill-down: compare exact methods pairwise to isolate framework vs solver cost.
Key question: if all RSib methods are equally slow, the bottleneck is the
shared framework (enumerate/rCall), NOT the solvers.
"""
import subprocess, random, tempfile, os, re
from collections import defaultdict

BINARY = os.path.join(os.path.dirname(__file__), "build", "release_build", "bk_algorithm")

def gen_and_write(n, p, seed, path):
    random.seed(seed)
    adj = defaultdict(set)
    edges = 0
    for u in range(n):
        for v in range(u+1, n):
            if random.random() < p:
                adj[u].add(v); adj[v].add(u); edges += 1
    with open(path, "w") as f:
        f.write(f"{n} {edges}\n")
        for v in range(n):
            nbrs = sorted(adj[v])
            f.write(str(v) + (" " + " ".join(map(str, nbrs)) if nbrs else "") + "\n")
    return n, edges

def run(gfile, mode, timeout=60):
    try:
        r = subprocess.run([BINARY, gfile, str(mode)],
                           capture_output=True, text=True, timeout=timeout)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        return None, None
    c = re.search(r"Total Maximal Cliques Found:\s*(\d+)|cliques=(\d+)", out)
    t = re.search(r"Time:\s*([\d.]+)\s*ms|time=([\d.]+)\s*ms", out)
    cliques = int(c.group(1) or c.group(2)) if c else None
    ms = float(t.group(1) or t.group(2)) if t else None
    return cliques, ms

REPS = 8

print("=" * 90)
print("DRILL-DOWN: Solver cost vs framework cost")
print("=" * 90)

# Test 1: How much does method choice matter at each graph size?
print("\n[A] Spread between fastest and slowest RSib method (ms) — shows solver cost")
print("    If spread is small → bottleneck is framework, not solver\n")
configs = [(20, 0.4), (35, 0.3), (50, 0.25), (75, 0.15), (100, 0.12)]
RSIB_MODES = {6: "BruteForce", 7: "Backtrack", 9: "Bitmask", 10: "MinHitSet"}

print(f"{'Config':<16} {'PivotBK':>10} {'BruteF':>10} {'Backtrack':>10} {'Bitmask':>10} {'MinHit':>10} {'Spread':>10} {'Spread%':>10}")
print("-" * 96)

with tempfile.TemporaryDirectory() as tmp:
    for n, p in configs:
        pivot_times, method_times = [], defaultdict(list)
        for rep in range(REPS):
            gfile = os.path.join(tmp, f"g.txt")
            gen_and_write(n, p, rep * 997 + n * 13, gfile)
            _, pt = run(gfile, 0)
            if pt: pivot_times.append(pt)
            for mode in RSIB_MODES:
                _, mt = run(gfile, mode)
                if mt: method_times[mode].append(mt)

        pavg = sum(pivot_times)/len(pivot_times) if pivot_times else 0
        avgs = {m: sum(method_times[m])/len(method_times[m]) for m in RSIB_MODES if method_times[m]}
        if avgs:
            mn, mx = min(avgs.values()), max(avgs.values())
            spread = mx - mn
            spread_pct = 100 * spread / mx if mx > 0 else 0
            row = f"n={n:<3} p={p:<5} {pavg:>10.3f}"
            for mode in RSIB_MODES:
                row += f" {avgs.get(mode, 0):>10.3f}"
            row += f" {spread:>10.3f} {spread_pct:>9.1f}%"
            print(row)

# Test 2: Solver-free estimate — what would RSib cost with NO solver?
# We compare RSib on a star graph (barely any cliques → collectCoveringCliques
# returns empty almost always → solver never runs) vs a dense graph.
print("\n\n[B] Solver vs framework isolation: sparse graph (solver rarely invoked)")
print("    vs dense graph (solver invoked often)")
print("    Framework cost ≈ sparse graph time; solver cost ≈ dense - sparse\n")

sparse_configs = [(50, 0.05), (100, 0.05), (150, 0.05)]
dense_configs  = [(50, 0.25), (100, 0.12), (150, 0.08)]

print(f"{'Config':<22} {'PivotBK':>10} {'BruteF':>10} {'Backtrack':>10} {'Bitmask':>10}")
print("-" * 65)

with tempfile.TemporaryDirectory() as tmp:
    for (ns, ps), (nd, pd) in zip(sparse_configs, dense_configs):
        for label, n, p in [("sparse", ns, ps), ("dense ", nd, pd)]:
            times = defaultdict(list)
            for rep in range(REPS):
                gfile = os.path.join(tmp, "g.txt")
                gen_and_write(n, p, rep * 997 + n, gfile)
                for mode in [0, 6, 7, 9]:
                    _, t = run(gfile, mode)
                    if t: times[mode].append(t)
            row = f"n={n} p={p} ({label})"
            row = f"{row:<22}"
            for mode in [0, 6, 7, 9]:
                avg = sum(times[mode])/len(times[mode]) if times[mode] else 0
                row += f" {avg:>10.3f}"
            print(row)
        print()

# Test 3: alreadyFound linear scan cost — compare graphs with few vs many cliques
print("\n[C] alreadyFound linear scan: few cliques vs many cliques at same n")
print("    More cliques = more overhead from find(allCliques.begin(),...)\n")

print(f"{'Config':<25} {'PivotBK':>10} {'BruteF':>10} {'#cliques':>10}")
print("-" * 60)

clique_configs = [
    (30, 0.2, "low density"),
    (30, 0.4, "mid density"),
    (30, 0.6, "high density"),
    (40, 0.2, "low density"),
    (40, 0.4, "mid density"),
]

with tempfile.TemporaryDirectory() as tmp:
    for n, p, dlabel in clique_configs:
        times_pivot, times_rsib, clique_counts = [], [], []
        for rep in range(REPS):
            gfile = os.path.join(tmp, "g.txt")
            gen_and_write(n, p, rep * 997 + n, gfile)
            c0, t0 = run(gfile, 0)
            c6, t6 = run(gfile, 6)
            if t0: times_pivot.append(t0)
            if t6: times_rsib.append(t6)
            if c0: clique_counts.append(c0)
        pavg = sum(times_pivot)/len(times_pivot) if times_pivot else 0
        ravg = sum(times_rsib)/len(times_rsib) if times_rsib else 0
        cavg = sum(clique_counts)/len(clique_counts) if clique_counts else 0
        label = f"n={n} p={p} {dlabel}"
        print(f"{label:<25} {pavg:>10.3f} {ravg:>10.3f} {cavg:>10.0f}")

print("\n" + "=" * 90)
print("CONCLUSIONS")
print("=" * 90)
print("""
Key questions to answer:
  [A] If spread between RSib methods is < 10% of total time → solver is NOT bottleneck
  [B] If sparse (solver rarely runs) ≈ dense (solver always runs) → framework IS bottleneck
  [C] If time grows with clique count at fixed n → alreadyFound scan is a bottleneck
""")
