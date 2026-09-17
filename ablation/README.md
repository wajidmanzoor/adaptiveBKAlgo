# Recursive-state ablation branch (historical `ablation-pxr-states`)

This branch compares implementation-native recursive search entries across
four exact maximal-clique systems:

- final CCRMCE-based reorder: the sum of FindOne and exhaustive CCRMCE
  recursive entries;
- HBBMC: `counter.vertex_recursive_calls` for its recursive `(R,P,X)` search;
- sparse adjacency-list Tomita: one state for every entry to
  `listAllMaximalCliquesAdjacencyListRecursive`;
- standalone CCRMCE: one state for every entry to `BKTomitaRecCCRV3`.

The branch name is retained for continuity, although the final reorder no
longer contains the old PXR implementation. These counters measure each
implementation's own recursive search state. They deliberately exclude input
conversion, ordering, preprocessing/root setup, reorder worklist and
hitting-set work, graph reduction, and terminal continuation expansion.
Because the algorithms decompose roots differently, the counts are useful as
implementation-level work measurements rather than identical abstract nodes.

All four systems use minimum clique size 3 and retain the actual vertex list
for every reported clique. Reorder uses budget 1000, the production pruning
profile, hit-set capacity 128, small-Q CCRMCE threshold 32, adaptive direct
threshold 256, and disabled ET1/ET2/ET3. HBBMC runs with
`--graph-reduction none --et 0`. Tomita uses its sparse adjacency-list
algorithm, and standalone CCRMCE uses CoreCliqueRemovalV3.

Run from the repository root:

```bash
./ablation/pxr_states/run.sh DATA_ROOT [RESULT_ROOT]
```

Useful controls are:

```text
PXR_REORDER_BUDGET=1000
PXR_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
PXR_TIMEOUT_SECONDS=3600
PXR_BUILD_JOBS=4
PXR_FAIL_ON_ERROR=0
```

Paired inputs must exist below `DATA_ROOT/adjacencylist/GROUP/GRAPH` and
`DATA_ROOT/edgelist/GROUP/GRAPH`. The runner converts the normalized edge list
to Tomita format before timing, records each system state count and ratio to
reorder, validates all four configurations and stored-list counts, and
requires clique-count agreement. Generated builds, converted inputs, logs, and
results are not tracked.
