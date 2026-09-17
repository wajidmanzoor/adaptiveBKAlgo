# Memory ablation

This runner compares peak resident memory for four exact maximal-clique
enumerators while every implementation retains the complete list of clique
vertex IDs. All systems report only maximal cliques of size at least 3.

The four systems are:

- `reorder`: the final CCRMCE-based reorder implementation, using the production
  pruning profile, threshold 32, and ET1/ET2/ET3 disabled;
- `hbbmc`: HBBMC++ with RMCE graph reduction and ET level 3;
- `tomita`: the sparse adjacency-list Tomita implementation compiled with
  `RETURN_CLIQUES_ONE_BY_ONE`;
- `ccrmce`: the standalone CoreCliqueRemovalV3 implementation.

## Run

From the repository root:

~~~bash
MEMORY_INTERVAL_US=1000 \
  ./ablation/memory/run.sh DATA_ROOT [RESULT_ROOT]
~~~

`DATA_ROOT` must contain matching `adjacencylist/GROUP/GRAPH` and
`edgelist/GROUP/GRAPH` files. Output defaults to
`ablation/memory/results/TIMESTAMP/`. Each mirrored group contains
`results.csv`, graph-specific logs, and traces. Reuse an explicit result
directory to resume.

The default `MEMORY_SYSTEMS=all` runs `reorder`, `hbbmc`, `tomita`, and
`ccrmce` sequentially. The legacy names `comparison` and `normal` are aliases
for the same four-system run. A comma-separated subset is also accepted, for
example `MEMORY_SYSTEMS=reorder,tomita`. Reorder is always run first when it
is selected so its clique count can serve as the reference.

`MEMORY_REORDER_BUDGET` changes the reorder seed-solver work budget from its
default of 1000; `unlimited` is accepted. `MEMORY_DATASETS` can select group
names, graph filenames, filename stems, or relative `GROUP/GRAPH` paths. An
empty value selects the complete corpus.

## Retained-clique contract

This study measures real clique storage, not a numeric counter approximation:

- reorder retains one explicit vertex vector per reported clique;
- HBBMC retains its direct-reduction and residual clique vectors;
- Tomita allocates an independent integer array for each clique and stores the
  arrays in its linked list;
- standalone CCRMCE stores every restored clique in
  `std::vector<std::vector<ui>>`.

Each executable prints both its enumerated and stored counts. The runner
requires `stored_cliques == cliques`, validates the intended algorithm
configuration, and compares HBBMC, Tomita, and CCRMCE counts with reorder. A
mismatch marks the campaign as failed when `MEMORY_FAIL_ON_ERROR=1`.

Tomita requires a symmetric comma-separated input. The runner converts the
paired normalized edge list before launching the sampler. That conversion is
recorded below `tomita_inputs/` and deliberately excluded from the measured
algorithm process.

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
`exec`, including graph parsing and retained clique storage; the sampler
process itself is excluded.

Each group `results.csv` summarizes count equality, trace length, sampled
peaks, observed `VmHWM`, and the final `wait4` peak. Full command, stdout, and
stderr files are kept below `GROUP/logs/GRAPH/`. Existing system/graph rows
with a corresponding trace are skipped when a result directory is resumed.
