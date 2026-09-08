# Pure seed-mask capacity ablation

This self-contained experiment compares fixed capacities `64`, `128`, `512`,
and `1024` with a dynamic mask. It builds only the physical Pure source copy
in this directory and never modifies or builds the production source in
`../../our/`.

All variants keep the normal Pure configuration fixed:

- optimized seed solver;
- budget 1000;
- ET1, ET2, and ET3 enabled;
- all seven pruning rules enabled;
- minimum clique size 3;
- small-Q full-PXR threshold 4;
- all cliques retained in the output list.

A fixed capacity `C` uses `C / 64` words and falls back to exact full PXR when
the reduced constraint count exceeds `C`. The dynamic variant allocates
exactly `ceil(reduced_constraints / 64)` words for each nontrivial seed-solver
call and therefore has no capacity fallback.

The common direct compact-arena path is held to one word (at most 64 raw
constraints) for every variant. Larger calls receive the same normalization
and unit-propagation pass before the chosen fixed or dynamic representation is
applied, so mask capacity is the only varying control.

Budget fallbacks and capacity fallbacks are reported separately. The fixed
budget of 1000 is intentional: this experiment selects the best capacity for
the production Pure configuration rather than measuring capacity under an
unlimited planning budget.

## Run

From the repository root:

~~~bash
./ablation/capacity/run.sh DATA_ROOT [RESULT_ROOT]
~~~

Pure inputs are discovered as `DATA_ROOT/adjacencylist/GROUP/GRAPH`. Output
defaults to `ablation/capacity/results/TIMESTAMP/`; every group receives
`results.csv`, `summary.csv`, and graph-specific logs. Reuse an explicit
result directory to resume.

The default campaign uses every discovered graph and three repetitions.

~~~text
CAPACITY_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
CAPACITY_REPETITIONS=3
CAPACITY_TIMEOUT_SECONDS=3600
CAPACITY_BUILD_JOBS=4
CAPACITY_FAIL_ON_ERROR=1
~~~

Each `results.csv` contains every raw run in that group. `summary.csv` reports
the geometric-mean wall-time ratio and speedup relative to capacity 128, mean
resource use, fallback totals, and exact-count agreement. A speedup above 1
favors the selected capacity over 128.

The runner records executable hashes, commands, stdout, stderr, and
`/usr/bin/time` resource measurements. Existing capacity/graph/repetition rows
are skipped when resuming a result directory.
