# Seed-mask capacity ablation branch

This branch contains only the seed-mask capacity study. Its physical source
copy is the final CCRMCE-based reorder and compares fixed capacities `64`,
`128`, `512`, and `1024` with a dynamic mask. Budget 1000, disabled legacy ET
terminals, the production pruning profile, adjacency hash threshold 256,
small-Q CCRMCE threshold 32, adaptive direct threshold 256, and minimum clique
size 3 remain fixed.

Run from the repository root:

```bash
./ablation/capacity/run.sh DATA_ROOT [RESULT_ROOT]
```

See [capacity/README.md](capacity/README.md) for representation semantics,
controls, CSV fields, exact-count validation, and resume behavior. Generated
builds, logs, and results are not tracked.
