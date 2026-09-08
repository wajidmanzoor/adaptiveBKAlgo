# Seed-mask capacity ablation branch

This branch contains only the seed-mask capacity study. Its physical Pure
source copy uses the latest optimized reorder implementation and compares
fixed capacities `64`, `128`, `512`, and `1024` with a dynamic mask. Budget
1000, ET1/ET2/ET3, all seven pruning rules, adjacency hash threshold 64,
small-Q full-PXR threshold 4, and minimum clique size 3 remain fixed.

Run from the repository root:

```bash
./ablation/capacity/run.sh DATA_ROOT [RESULT_ROOT]
```

See [capacity/README.md](capacity/README.md) for representation semantics,
controls, CSV fields, exact-count validation, and resume behavior. Generated
builds, logs, and results are not tracked.
