#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path


def erdos_renyi_edges(n: int, p: float, rng: random.Random):
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    log_q = math.log1p(-p)
    v = 1
    w = -1
    while v < n:
        r = rng.random()
        w += 1 + int(math.log1p(-r) / log_q)
        while w >= v and v < n:
            w -= v
            v += 1
        if v < n:
            yield v, w


def write_graph(path: Path, n: int, p: float, seed: int) -> tuple[int, int]:
    rng = random.Random(seed)
    adj = [[] for _ in range(n)]
    m = 0
    for u, v in erdos_renyi_edges(n, p, rng):
        adj[u].append(v)
        adj[v].append(u)
        m += 1

    for nbrs in adj:
        nbrs.sort()

    with path.open("w", encoding="utf-8") as f:
        f.write(f"{n} {m}\n")
        for u, nbrs in enumerate(adj):
            if nbrs:
                f.write(f"{u} {' '.join(map(str, nbrs))}\n")
            else:
                f.write(f"{u}\n")
    return n, m


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate random Erdos-Renyi graph datasets in the adjacency-list format used by bk_algorithm."
    )
    parser.add_argument(
        "--out-dir",
        default="data",
        help="Output directory for generated graph files (default: data)",
    )
    parser.add_argument(
        "--sparse-name",
        default="generated_random_sparse_n12000_p0.0015_seed20260525.txt",
    )
    parser.add_argument("--sparse-n", type=int, default=12000)
    parser.add_argument("--sparse-p", type=float, default=0.0015)
    parser.add_argument("--sparse-seed", type=int, default=20260525)
    parser.add_argument(
        "--dense-name",
        default="generated_random_dense_n6500_p0.008_seed20260525.txt",
    )
    parser.add_argument("--dense-n", type=int, default=6500)
    parser.add_argument("--dense-p", type=float, default=0.008)
    parser.add_argument("--dense-seed", type=int, default=20260525)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        (args.sparse_name, args.sparse_n, args.sparse_p, args.sparse_seed, "sparse"),
        (args.dense_name, args.dense_n, args.dense_p, args.dense_seed, "dense"),
    ]

    for name, n, p, seed, label in configs:
        path = out_dir / name
        n_written, m_written = write_graph(path, n, p, seed)
        print(
            f"{label}: wrote {path} with n={n_written}, m={m_written}, p={p}, seed={seed}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
