# Memory ablation

This runner compares normal Pure against HBBMC++ while both programs retain
their complete maximal-clique lists.

## Run

From the repository root:

~~~bash
MEMORY_INTERVAL_US=1000 \
  ./ablation/memory/run.sh DATA_ROOT [RESULT_ROOT]
~~~

`DATA_ROOT` must contain matching
`adjacencylist/GROUP/GRAPH` and `edgelist/GROUP/GRAPH` files. Output defaults
to `ablation/memory/results/TIMESTAMP/`. Each mirrored group contains
`results.csv`, graph-specific logs, and traces. Reuse an explicit result
directory to resume.

The default `comparison` profile runs exactly these systems sequentially:

- `pure`: optimized Pure, ET1/ET2/ET3 enabled, budget 1000, small-Q
  full-PXR threshold 4, and all seven pruning rules enabled;
- `hbbmc_plus_plus`: RMCE graph reduction and ET level 3.

`MEMORY_PURE_BUDGET` changes the Pure budget when a separate budget comparison
is intentional. The runner reads Pure's printed configuration and rejects a
run if ET, the selected budget, or any pruning rule is not enabled as expected.
It likewise rejects the HBBMC++ run unless RMCE, ET 3, and the
`hbbmc_plus_plus=true` configuration are confirmed.

Set `MEMORY_SYSTEMS=all` to add unreduced ET 1/2/3 and reduced ET 1/2, covering
all eight HBBMC runtime combinations. Individual comma-separated names are
also accepted:

~~~text
pure,hbbmc,hbbmc_et1,hbbmc_et2,hbbmc_et3,hbbmc_plus,
hbbmc_gr_et1,hbbmc_gr_et2,hbbmc_plus_plus
~~~

`MEMORY_DATASETS` can select group names, graph filenames, filename stems, or
relative `GROUP/GRAPH` paths. An empty value selects the complete corpus.

## Trace format

Each run writes `GROUP/traces/GRAPH/SYSTEM.txt`:

~~~text
# interval_us=1000
# columns=elapsed_us rss_kb vm_hwm_kb
0       352     352
1058    4216    4216
...
# samples=...
# peak_sampled_rss_kb=...
# peak_observed_hwm_kb=...
# wait4_peak_rss_kb=...
# timed_out=0
~~~

Sampling uses `/proc/PID/status` and therefore requires Linux. The sampler
uses monotonic deadlines, but OS scheduling cannot guarantee execution at an
exact microsecond boundary. `elapsed_us` records when each sample actually
occurred. The trace covers the direct algorithm process from successful
`exec`, including graph input parsing and retained output storage; the sampler
process itself is excluded.

Each group `results.csv` summarizes count equality, trace length, sampled
peaks, the observed `VmHWM`, and the final `wait4` peak. Full command, stdout,
and stderr files are kept below `GROUP/logs/GRAPH/`. Existing system/graph rows
with a corresponding trace are skipped when a result directory is resumed.
