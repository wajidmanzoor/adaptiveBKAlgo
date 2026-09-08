# Keep-one pruning ablation branch

This branch contains only the seven-rule keep-one pruning study. Its physical
Pure source copy uses the latest optimized reorder implementation. The runner
builds an all-rules reference plus variants retaining exactly one of
normalization, subsumption, unit propagation, usefulness, antichain,
fail-first, or zero-coverage pruning while disabling the other six.

Budget 1000, ET1/ET2/ET3, hit-set capacity 128, adjacency hash threshold 64,
small-Q full-PXR threshold 4, and minimum clique size 3 remain fixed.

Run from the repository root:

```bash
./ablation/pruning_keep_one/run.sh DATA_ROOT [RESULT_ROOT]
```

See [pruning_keep_one/README.md](pruning_keep_one/README.md) for controls,
validation, CSV output, and resume behavior. Generated builds, logs, and
results are not tracked.
