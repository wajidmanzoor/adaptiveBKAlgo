# Keep-one-rule pruning ablation

This complementary study measures how Pure behaves when exactly one of its
seven pruning rules remains enabled. It is separate from the leave-one-out
study in `../pruning/`.

## Run

~~~bash
./run.sh DATA_ROOT [RESULT_ROOT]
~~~

The runner builds an all-rules reference followed by seven variants:
`only_normalization`, `only_subsumption`, `only_unit`,
`only_usefulness`, `only_antichain`, `only_fail_first`, and
`only_zero_coverage`. In every `only_*` build, the named rule is enabled and
the other six are disabled.

Every completed run validates its compiled seven-rule mask, ET configuration,
budget, fixed small-Q threshold of 4, numeric output, and exact clique count
against the all-rules reference. Runs are sequential and resumable.

## Controls

~~~text
KEEP_ONE_BUDGET=1000
KEEP_ONE_DATASETS=GROUP,GRAPH,OR_RELATIVE_PATH
KEEP_ONE_TIMEOUT_SECONDS=3600
KEEP_ONE_BUILD_JOBS=4
KEEP_ONE_FAIL_ON_ERROR=0
~~~

By default, results are written below `results/TIMESTAMP/`. Supplying an
explicit result root resumes that campaign. Each input group receives a
`results.csv`; commands, stdout, stderr, timing/resource records, build logs,
binaries, hashes, and environment metadata are retained.
