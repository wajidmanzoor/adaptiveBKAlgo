# Memory ablation branch

This branch isolates the memory study around the final CCRMCE-based reorder
implementation. It compares four exact systems: final reorder, HBBMC++, the
sparse adjacency-list Tomita implementation, and standalone CCRMCE. Each
system retains the actual vertex list for every maximal clique of size at
least 3, and the runner rejects stored-count or cross-system count mismatches.

Run from the repository root:

```bash
MEMORY_INTERVAL_US=1000 \
  ./ablation/memory/run.sh DATA_ROOT [RESULT_ROOT]
```

Paired inputs must exist below `DATA_ROOT/adjacencylist/GROUP/GRAPH` and
`DATA_ROOT/edgelist/GROUP/GRAPH`. See
[memory/README.md](memory/README.md) for system selection, retained-list
semantics, trace format, validation, and resume behavior. Generated builds,
converted inputs, traces, logs, and results are not tracked.
