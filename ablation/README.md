# Leave-one-out pruning ablation branch

This branch contains only the seven-rule leave-one-out pruning study. Its
physical Pure source copy uses the latest optimized reorder implementation.
The runner builds an all-rules reference plus variants disabling exactly one
of normalization, subsumption, unit propagation, usefulness, antichain,
fail-first, or zero-coverage pruning.

Budget 1000, ET1/ET2/ET3, hit-set capacity 128, adjacency hash threshold 64,
small-Q full-PXR threshold 4, and minimum clique size 3 remain fixed.

Run from the repository root:

```bash
./ablation/pruning/run.sh DATA_ROOT [RESULT_ROOT]
```

Controls:

```text
PRUNING_BUDGET=1000
PRUNING_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
PRUNING_TIMEOUT_SECONDS=3600
PRUNING_BUILD_JOBS=4
PRUNING_FAIL_ON_ERROR=0
```

Every completed build validates its compiled rule mask and exact clique count
against the all-rules reference. Generated builds, logs, and results are not
tracked.
