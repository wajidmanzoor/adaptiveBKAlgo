# Recursive PXR-state ablation branch

This branch contains only the recursive PXR-state comparison. It includes
physical Pure and HBBMC source copies. Pure uses the latest optimized reorder
core with ET1/ET2/ET3 disabled at compile time; HBBMC runs with
`--graph-reduction none --et 0`. Both retain all reported cliques.

Pure counts recursive `(R,P,X)` entries in FindOne and full-PXR fallback.
HBBMC uses `counter.vertex_recursive_calls`. Ordering, worklist entries,
seed-solver work, graph reduction, and terminal continuation work are excluded.
Pure budget 1000, all seven pruning rules, hit-set capacity 128, adjacency hash
threshold 64, small-Q full-PXR threshold 4, and minimum clique size 3 remain
fixed.

Run from the repository root:

```bash
./ablation/pxr_states/run.sh DATA_ROOT [RESULT_ROOT]
```

Paired inputs must exist below `DATA_ROOT/adjacencylist/GROUP/GRAPH` and
`DATA_ROOT/edgelist/GROUP/GRAPH`. The runner validates both configurations,
state-count invariants, stored clique counts, and cross-implementation clique
count equality. Generated builds, logs, and results are not tracked.
