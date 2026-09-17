# Recursive-state ablation branch (historical `ablation-pxr-states`)

This branch compares recursive search states in the final CCRMCE-based reorder
against HBBMC. The branch name is retained for continuity, but the final
reorder no longer contains PXR: its measured state is one recursive CCRMCE
entry in FindOne or exhaustive enumeration. HBBMC continues to use
`counter.vertex_recursive_calls` for its recursive `(R,P,X)` entries.

Both physical source copies retain every reported clique. Reorder uses budget
1000, the production pruning profile, hit-set capacity 128, adjacency hash
threshold 256, small-Q CCRMCE threshold 32, adaptive direct threshold 256,
and disabled ET1/ET2/ET3. HBBMC runs with
`--graph-reduction none --et 0`. Ordering, worklist entries, seed-solver work,
graph reduction, and terminal continuation work are excluded.

Run from the repository root:

```bash
./ablation/pxr_states/run.sh DATA_ROOT [RESULT_ROOT]
```

Paired inputs must exist below `DATA_ROOT/adjacencylist/GROUP/GRAPH` and
`DATA_ROOT/edgelist/GROUP/GRAPH`. The runner validates both configurations,
state-count invariants, stored clique counts, and cross-implementation clique
count equality. Generated builds, logs, and results are not tracked.
