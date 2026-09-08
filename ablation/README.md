# Budget ablation branch

This branch contains only the budget study. Its physical Pure source copy uses
the latest optimized reorder implementation while varying the per-seed solver
work budget. ET1/ET2/ET3, all seven pruning rules, hit-set capacity 128,
adjacency hash threshold 64, small-Q full-PXR threshold 4, and minimum clique
size 3 remain fixed.

Run from the repository root:

```bash
./ablation/budgets/run.sh DATA_ROOT [RESULT_ROOT]
```

Pure inputs are discovered below `DATA_ROOT/adjacencylist/GROUP/GRAPH`.
Results default to `ablation/budgets/results/TIMESTAMP/` and are resumable when
an explicit result directory is reused.

Controls:

```text
BUDGET_VALUES=0,100,1000,10000,100000,unlimited
BUDGET_REFERENCE=10000
BUDGET_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
BUDGET_TIMEOUT_SECONDS=3600
BUDGET_BUILD_JOBS=4
BUDGET_FAIL_ON_ERROR=0
```

Every completed variant validates its runtime configuration and exact clique
count against the reference budget. Generated builds, logs, and results are
not tracked.
