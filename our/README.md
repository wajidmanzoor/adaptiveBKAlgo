# CCRMCE inside AdaptiveBK/ReorderSib

This directory is an experimental copy of the GitHub-main
AdaptiveBK/ReorderSib implementation. Its reorder, covering-clique, exact
sibling-seed, work-budget, and worklist logic remain intact. Every branch
search lane formerly handled by PXR now uses Core Clique Removal MCE:

- `findOnePure` uses stop-after-one CCRMCE;
- an over-budget seed solver uses exhaustive CCRMCE; and
- every `|Q| <= 4` terminal uses a one-word exhaustive CCRMCE kernel.

For a formal branch `B=(M,Q)`, minimum-degree peeling partitions `Q` into
a core clique `C` and residual candidates `P`. Recursion represents
`C/P/X` over `Q` as bitsets, keeps external `X` compact, pivots over
`C union P union X`, and uses the CCRMCE `T1 union T2` branch set. The
find-one lane filters external `X` lazily; the exhaustive lane materializes
its `X`-to-`Q` rows once. Local adjacency construction chooses between CSR
scans and pairwise lookups according to estimated work.

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
It runs minimum clique sizes 1 and 3; budget zero exercises exhaustive
CCRMCE, while budget 1000 exercises normal find-one and sibling-seed routing.

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

Final optimized results are in
`results/optimized_ccr_pxr_20260913/`. The earlier
`results/ccr_pxr_20260913/` directory is retained as the documented
correctness-first vector implementation; it is not representative of the
current bitset engine.
