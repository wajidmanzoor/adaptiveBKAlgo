# CCRMCE inside AdaptiveBK/ReorderSib

This directory is an experimental copy of the GitHub-main
AdaptiveBK/ReorderSib implementation. Its reorder, covering-clique, exact
sibling-seed, work-budget, and worklist logic remain intact. Every branch
search lane formerly handled by PXR now uses Core Clique Removal MCE:

- `findOnePure` uses stop-after-one CCRMCE;
- an over-budget seed solver uses exhaustive CCRMCE; and
- worklist branches with `|Q| <= 32` route directly to exhaustive CCRMCE.

The exhaustive lane uses a materialized scalar one-word kernel through 64
candidates and a fixed two-word kernel for the retained 68--128 range. On
productive roots, adaptive direct routing extends to 256 candidates. Wider
branches still use the unchanged covering-clique and exact sibling-seed path.

For a formal branch `B=(M,Q)`, minimum-degree peeling partitions `Q` into
a core clique `C` and residual candidates `P`. Recursion represents
`C/P/X` over `Q` as bitsets, keeps external `X` compact, pivots over
`C union P union X`, and uses the CCRMCE `T1 union T2` branch set. The
find-one lane filters external `X` lazily; the exhaustive lane materializes
its `X`-to-`Q` rows once. Local adjacency construction chooses between CSR
scans and pairwise lookups according to estimated work. Global adjacency
queries use sorted CSR for low-degree rows and a compact flat hash for rows
of degree at least 256.

ET1, ET2, and ET3 are disabled. The retained FastPlex3 code is not called by
these CCRMCE lanes.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The correctness test checks every undirected graph through six vertices plus
600 fixed-seed random graphs on 7--12 vertices against a brute-force oracle.
It runs minimum clique sizes 1 and 3 at budgets 0 and 1000. A separate
327-vertex graph with four analytically known maximal cliques forces the
positive-budget exact sibling-seed route and crosses the flat-hash threshold.

## Run

```bash
./build/adaptive_bk PATH/TO/GRAPH --budget 1000 --min-clique-size 3
```

The output includes `reorder.ccr.*` counters for core extraction and both
CCRMCE lanes.

## Benchmark

`scripts/run_external_comparison.sh` records status, exact clique count,
internal enumeration time, wall time, and peak RSS for one or more binaries.
`scripts/summarize_comparison.py` validates completed counts and produces
paired medians and speedups.

The final full-corpus run is in
`results/reorder_ccrmce_optimized_20260914.csv`. The earlier
optimized sequence is in `results/optimized_ccr_pxr_20260913/`, while the
`results/ccr_pxr_20260913/` directory is retained as the documented
correctness-first vector implementation; it is not representative of the
current bitset engine.
