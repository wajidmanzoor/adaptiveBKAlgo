# Budget ablation branch

This branch contains only the budget study. Its physical source copy is the
final CCRMCE-based reorder implementation from `latestUpdate`, augmented only
with seed-solver fallback counters while varying the per-call work budget.
ET1/ET2/ET3 remain disabled, hit-set capacity remains 128, the small-branch
CCRMCE threshold remains 32, the adaptive direct threshold remains 256, and
minimum clique size remains 3. The production pruning profile is preserved:
subsumption is disabled and the other six rules are enabled.

Run from the repository root:

```bash
./ablation/budgets/run.sh DATA_ROOT [RESULT_ROOT]
```

Reorder inputs are discovered below `DATA_ROOT/adjacencylist/GROUP/GRAPH`.
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

Every completed variant validates the final engine configuration, verifies
that no capacity fallback occurred, and checks its exact clique count against
the reference budget. Generated builds, logs, and results are not tracked.
