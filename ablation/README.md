# Leave-one-out pruning ablation branch

This branch contains only the seven-rule leave-one-out pruning study. Its
physical source copy is the final CCRMCE-based reorder from `latestUpdate`,
augmented with compile-time rule switches and measurement counters. The runner
builds an all-rules reference plus variants disabling exactly one of
normalization, subsumption, unit propagation, usefulness, antichain,
fail-first, or zero-coverage pruning.

The all-rules reference deliberately enables subsumption so every rule can be
removed symmetrically; production `latestUpdate` otherwise leaves subsumption
disabled. Budget 1000, disabled ET1/ET2/ET3, hit-set capacity 128, adjacency
hash threshold 256, small-Q CCRMCE threshold 32, adaptive direct threshold
256, and minimum clique size 3 remain fixed.

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
