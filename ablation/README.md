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
for every reported clique. Reorder uses an unlimited solver budget, the
production pruning profile, hit-set capacity 128, small-Q CCRMCE threshold 32,
adaptive direct threshold 256, and disabled ET1/ET2/ET3. HBBMC runs with
`--graph-reduction none --et 0`. Tomita uses its sparse adjacency-list
algorithm, and standalone CCRMCE uses CoreCliqueRemovalV3.

The branch-local binary in `ablation/pxr_states/pure/` also defaults to
`unlimited` when `--budget` is omitted. The artifact runner passes
`--budget unlimited` explicitly, records it in `environment.txt`, and rejects
a run whose output reports any other value.

Run from the repository root:

```bash
./ablation/pxr_states/run.sh ADJACENCY_DIR EDGE_DIR [RESULT_ROOT]
```

Useful controls are:

```text
PXR_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
PXR_TIMEOUT_SECONDS=1800
PXR_BUILD_JOBS=4
PXR_FAIL_ON_ERROR=0
```

Paired inputs must have the same relative `GROUP/GRAPH` path below
`ADJACENCY_DIR` and `EDGE_DIR`. The runner converts the normalized edge list
to Tomita format before timing, records each system state count and ratio to
reorder, validates all four configurations and stored-list counts, and
requires clique-count agreement. A failed conversion or algorithm run is
recorded and later graphs continue. Generated builds, converted inputs, logs,
and results are not tracked.
