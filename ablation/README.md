# Memory ablation branch

This branch contains only the memory study. It includes physical copies of the
latest optimized reorder implementation and the HBBMC++ comparison, plus the
Linux RSS sampler. The default comparison runs Pure with budget 1000,
ET1/ET2/ET3, all seven pruning rules, hit-set capacity 128, adjacency hash
threshold 64, small-Q full-PXR threshold 4, and minimum clique size 3 against
HBBMC++ with RMCE reduction and ET level 3.

Run from the repository root:

```bash
MEMORY_INTERVAL_US=1000 \
  ./ablation/memory/run.sh DATA_ROOT [RESULT_ROOT]
```

Paired inputs must exist below `DATA_ROOT/adjacencylist/GROUP/GRAPH` and
`DATA_ROOT/edgelist/GROUP/GRAPH`. See
[memory/README.md](memory/README.md) for selectable systems, trace format,
validation, and resume behavior. Generated builds, traces, logs, and results
are not tracked.
