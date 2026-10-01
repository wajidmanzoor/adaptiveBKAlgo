# VLDB Artifact: AdaptiveBK / ReorderSib

This repository is the artifact for an exact maximal-clique enumeration system
based on the optimized ReorderSib worklist algorithm. The `main` branch contains
the final implementation and the primary four-system wall-time experiment.
Instrumentation-heavy ablations are isolated on separate branches so they
cannot change the production binary used by the main experiment.

The final search combines ascending degeneracy reordering, CSR adjacency,
exact sibling generation through a 128-constraint hitting-set solver, a bounded
seed solver with exact CCRMCE fallback, specialized small-branch enumeration,
compact clique storage, and exact duplicate detection. All paper experiments
report maximal cliques of size at least three.

## Artifact map

| Evaluation target | Branch | Entry point | Primary output |
| --- | --- | --- | --- |
| Four-system end-to-end wall time | [`main`](https://github.com/wajidmanzoor/adaptiveBKAlgo/tree/main) | `scripts/run_benchmark.sh` | `runs.csv`, `comparison.csv` |
| Seed-solver work budget | [`ablation-budgets`](https://github.com/wajidmanzoor/adaptiveBKAlgo/tree/ablation-budgets) | `ablation/budgets/run.sh` | per-group `results.csv` |
| Seed-mask capacity | [`ablation-capacity`](https://github.com/wajidmanzoor/adaptiveBKAlgo/tree/ablation-capacity) | `ablation/capacity/run.sh` | per-group `results.csv`, `summary.csv` |
| Retained-clique memory | [`ablation-memory`](https://github.com/wajidmanzoor/adaptiveBKAlgo/tree/ablation-memory) | `ablation/memory/run.sh` | per-group `results.csv`, RSS traces |
| Recursive search states | [`ablation-pxr-states`](https://github.com/wajidmanzoor/adaptiveBKAlgo/tree/ablation-pxr-states) | `ablation/pxr_states/run.sh` | per-group `results.csv` |

Start with the quick validation below. The complete corpus includes graphs with
millions of vertices and graphs that produce very large clique collections, so
a full campaign is intended for a server rather than a laptop.

## 1. Requirements

The supported evaluation environment is 64-bit Linux. No GPU is required.

- CMake 3.16 or newer;
- GCC or Clang with C99 and C++17 support;
- Bash;
- GNU `awk`, `sort`, `find`, `timeout`, `date`, `sha256sum`, and
  `/usr/bin/time`;
- enough RAM and storage for the selected graph and its retained maximal
  cliques.

The memory ablation additionally reads `/proc/PID/status`, so it must run on
Linux. Release builds use `-O3` or CMake Release optimization and
`-march=native`; compare variants on the same machine. Every campaign is
sequential and launches each algorithm with `OMP_NUM_THREADS=1`.

For Ubuntu 24.04 or a similar Debian-based system, the core packages are:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake time gawk coreutils
```

## 2. Obtain and validate the artifact

```bash
git clone https://github.com/wajidmanzoor/adaptiveBKAlgo.git
cd adaptiveBKAlgo
git switch main

cmake -S our -B our/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON
cmake --build our/build --parallel 4
ctest --test-dir our/build --output-on-failure
```

The executable is `our/build/adaptive_bk`. The correctness test compares every
undirected graph through six vertices, 600 fixed-seed random graphs on 7--12
vertices, and one structured threshold-crossing graph with a brute-force
oracle.

Run one bundled graph directly:

```bash
./our/build/adaptive_bk data/real/as-733.txt \
  --budget 1000 \
  --min-clique-size 3
```

A successful run prints key/value records including:

```text
reorder.cliques=...
reorder.stored_cliques=...
reorder.minimum_clique_size=3
reorder.budget=1000
reorder.runtime_ms=...
```

`reorder.cliques` and `reorder.stored_cliques` must agree. Use
`--min-clique-size 1` only when conventional maximal-clique enumeration must
also include maximal singletons and two-vertex cliques. `--print-cliques`
prints the canonical original vertex IDs.

## 3. Graph files and provenance

### Bundled processed files

The following repository paths are links to ready-to-use Reorder-format graph
files:

- [`data/real/`](data/real/) contains the processed real-world corpus that is
  small enough to redistribute, plus the five processed random graphs;
- [`data/data_rand200/`](data/data_rand200/) contains 200 small random graphs
  used during development and validation;
- [`data/`](data/) contains additional synthetic and regression inputs;
- [`compare/HBBMC/edgeListCleaned/`](compare/HBBMC/edgeListCleaned/) contains
  small normalized edge-list examples.

The five processed random inputs used by the artifact are directly available
here:

| File | Vertices | Undirected edges |
| --- | ---: | ---: |
| [`randomGraph1.txt`](data/real/randomGraph1.txt) | 500 | 6,276 |
| [`randomGraph2.txt`](data/real/randomGraph2.txt) | 500 | 62,620 |
| [`randomGraph3.txt`](data/real/randomGraph3.txt) | 2,500 | 31,256 |
| [`randomGraph4.txt`](data/real/randomGraph4.txt) | 2,500 | 312,386 |
| [`randomGraph5.txt`](data/real/randomGraph5.txt) | 10,000 | 124,590 |

Large third-party graphs are not duplicated in Git. Obtain them from their
original collections and retain the source license and citation:

- [Stanford Large Network Dataset Collection (SNAP)](https://snap.stanford.edu/data/)
  for the AS, CA-HepTh, DBLP, Enron, Orkut, Patents, Pokec, Skitter, Stanford,
  and YouTube graphs;
- [Network Data Repository](https://networkrepository.com/networks.php) for
  additional social, web, communication, and biological networks;
- [SuiteSparse Matrix Collection](https://sparse.tamu.edu/) for the GSE,
  GR-QC, NASA, ship, road, and other sparse-matrix graphs;
- [GTgraph](https://www.cse.psu.edu/~kxm85/software/GTgraph/) for synthetic
  R-MAT graphs.

Names ending in `_q` or `_sort_q` are processed local names. Search the source
collection using the base name, then apply the normalization below. Some source
collections are directed, weighted, one-based, contain self-loops, or contain
parallel edges; the artifact always converts them to a simple undirected graph.

### Required paired input formats

The main, memory, and recursive-state runners use two directory trees with
identical relative paths:

```text
ADJACENCY_DIR/GROUP/GRAPH
EDGE_DIR/GROUP/GRAPH
```

The adjacency file begins with `n m` and then has exactly `n` rows in vertex
order. Each row starts with its vertex ID followed by its sorted neighbors:

```text
4 4
0 1 3
1 0 2
2 1 3
3 0 2
```

The paired edge file has no header and contains each undirected edge exactly
once with `u < v`:

```text
0 1
0 3
1 2
2 3
```

IDs must be integers in `0..n-1`. Isolates still require an adjacency row.
Self-loops are removed, direction is ignored, and duplicate or reciprocal
edges are collapsed.

### Process a raw edge list

[`scripts/prepare_graph.sh`](scripts/prepare_graph.sh) performs the exact
normalization with GNU external sorting, so it does not need to hold all edges
in a Python set. The raw input must have integer endpoints in its first two
whitespace-separated fields; later fields such as weights are ignored.

For a zero-based graph whose largest ID determines `n`:

```bash
./scripts/prepare_graph.sh \
  raw/example.edges \
  prepared/adjacency/real/example.graph \
  prepared/edges/real/example.graph
```

For a one-based graph with 10,000 vertices, including possible isolates:

```bash
GRAPH_VERTEX_BASE=1 ./scripts/prepare_graph.sh \
  raw/example.edges \
  prepared/adjacency/real/example.graph \
  prepared/edges/real/example.graph \
  10000
```

The helper rejects noninteger IDs and existing outputs. Set `GRAPH_FORCE=1`
only when replacement is intentional. If a source uses sparse or textual IDs,
first create and retain a deterministic mapping to contiguous integer IDs.
Always pass the source vertex count when isolated vertices must be preserved.

To derive the paired edge file from an already processed bundled adjacency
file, emit only the upper triangle and verify its line count against the
header:

```bash
mkdir -p prepared/adjacency/smoke prepared/edges/smoke
cp data/real/as-733.txt prepared/adjacency/smoke/as-733.txt
awk 'NR > 1 { u=$1; for (i=2; i<=NF; ++i) if (u < $i) print u, $i }' \
  data/real/as-733.txt > prepared/edges/smoke/as-733.txt

expected_edges=$(awk 'NR == 1 { print $2 }' data/real/as-733.txt)
actual_edges=$(wc -l < prepared/edges/smoke/as-733.txt)
test "$expected_edges" -eq "$actual_edges"
```

### Generate the random graphs with GTgraph R-MAT

The synthetic random graphs for this artifact use the R-MAT generator in
[GTgraph](https://www.cse.psu.edu/~kxm85/software/GTgraph/), not the uniform
`GTgraph-random` generator. The
[GTgraph generator guide](https://www.cse.psu.edu/~kxm85/software/GTgraph/gen.pdf)
documents `GTgraph-rmat`, its DIMACS output, and its configuration file.

After downloading and building GTgraph, generate an instance with the recorded
vertex and edge parameters:

```bash
./GTgraph-rmat -n N -m M -o raw/rmat-N-M.gr
```

GTgraph's documented default R-MAT quadrant probabilities are
`a=0.45`, `b=0.15`, `c=0.15`, and `d=0.25`. If a custom configuration is used,
pass it with `-c CONFIG` and archive the configuration with the results. Also
record the generator version, requested `N` and `M`, and seed. The final simple
undirected edge count can be smaller than requested `M` after self-loops,
parallel edges, and reciprocal edges are collapsed.

GTgraph writes one-based DIMACS arc records. Extract their endpoints, read the
declared vertex count from the `p` record, and create both artifact formats:

```bash
num_vertices=$(awk '$1 == "p" { print $3; exit }' raw/rmat-N-M.gr)
awk '$1 == "a" { print $2, $3 }' raw/rmat-N-M.gr \
  > raw/rmat-N-M.edges

GRAPH_VERTEX_BASE=1 ./scripts/prepare_graph.sh \
  raw/rmat-N-M.edges \
  prepared/adjacency/random/rmat-N-M.graph \
  prepared/edges/random/rmat-N-M.graph \
  "$num_vertices"
```

Use the bundled `randomGraph1.txt`--`randomGraph5.txt` files when exact reuse of
the processed artifact inputs is required. Regeneration is useful for extending
the study; byte-for-byte regeneration additionally depends on the archived
GTgraph seed and configuration.

## 4. Main experiment: four-system wall time

The primary experiment compares:

1. final Reorder/AdaptiveBK with budget 1,000;
2. the independent paper-faithful HBBMC++ implementation with RMCE reduction
   and ET level 3;
3. the sparse adjacency-list Tomita implementation; and
4. standalone CCRMCE.

The measured `wall_time_us` encloses each complete native process invocation:
startup, native-format parsing, preprocessing, enumeration, retained output
work, and teardown. Tomita's one-time format conversion is logged separately
and excluded from that interval. All four systems retain their reported clique
output, and the runner checks exact clique-count agreement.

### Fast artifact check

Prepare `as-733` as shown above, then run only that graph:

```bash
COMPARE_DATASETS=as-733 \
COMPARE_TIMEOUT_SECONDS=300 \
COMPARE_FAIL_ON_ERROR=1 \
  bash scripts/run_benchmark.sh \
    prepared/adjacency \
    prepared/edges \
    results/walltime-smoke
```

The script configures and builds all four implementations. Success ends with a
`CAMPAIGN_COMPLETE` line. In
`results/walltime-smoke/comparison.csv`, the `all_counts_match` value should be
`yes`.

### Full campaign

Place every selected graph at the same relative path in the two input roots,
then run:

```bash
COMPARE_TIMEOUT_SECONDS=1800 \
COMPARE_BUILD_JOBS=4 \
COMPARE_FAIL_ON_ERROR=1 \
  bash scripts/run_benchmark.sh \
    /absolute/path/to/adjacency \
    /absolute/path/to/edges \
    /absolute/path/to/results/walltime
```

Useful controls are:

```text
COMPARE_TIMEOUT_SECONDS=1800
COMPARE_BUILD_JOBS=4
COMPARE_DATASETS=GROUP,GRAPH,STEM,OR_RELATIVE_PATH
COMPARE_FAIL_ON_ERROR=0
```

An empty `COMPARE_DATASETS` selects the complete corpus. A timeout or failure
is recorded and later programs and graphs continue. Reusing the same result
directory resumes the campaign: an existing variant/graph row is not rerun.

The result directory contains:

```text
environment.txt          exact configuration, binary SHA-256 values, uname
runs.csv                 one row per graph and implementation
comparison.csv           aligned counts and wall times for all four systems
build_logs/              CMake configure and build logs
logs/GROUP/              command, stdout, stderr, and resource files
tomita_inputs/GROUP/     generated Tomita input and conversion logs
```

Interpret `wall_time_us` only for rows whose status is `completed`. Do not turn
timeouts into numeric runtimes. `comparison.csv` uses these wall-time columns:
Reorder 7, HBBMC 10, Tomita 13, and CCRMCE 16. For a quick per-graph speedup
table relative to Reorder:

```bash
awk -F, '
  BEGIN { print "graph,hbbmc_over_reorder,tomita_over_reorder,ccrmce_over_reorder" }
  NR > 1 && $5 == "completed" && $8 == "completed" &&
            $11 == "completed" && $14 == "completed" && $17 == "yes" {
    printf "%s,%.6f,%.6f,%.6f\n", $2, $10/$7, $13/$7, $16/$7
  }
' results/walltime/comparison.csv
```

A ratio above 1 means Reorder is faster than that comparison system. Use
geometric means for multiplicative speedups and keep the timeout set visible in
any reported aggregate.

## 5. Ablation experiments

Each ablation must be run from its named branch. Do not copy an ablation
directory onto `main`: the branch also contains matching common helpers and an
isolated source copy. Prepare the paired graph trees before switching branches,
prefer absolute input/result paths, and commit no generated results.

Fetch the branches once:

```bash
git fetch origin \
  ablation-budgets \
  ablation-capacity \
  ablation-memory \
  ablation-pxr-states
```

### 5.1 Work-budget ablation

```bash
git switch ablation-budgets
git pull --ff-only origin ablation-budgets

BUDGET_VALUES=0,100,1000,10000,100000,unlimited \
BUDGET_REFERENCE=10000 \
BUDGET_TIMEOUT_SECONDS=1800 \
BUDGET_FAIL_ON_ERROR=1 \
  ./ablation/budgets/run.sh \
    /absolute/path/to/adjacency \
    /absolute/path/to/results/budgets
```

This varies only the per-call seed-solver compatibility-work budget. Hit-set
capacity remains 128, minimum clique size remains 3, and the production
pruning and CCRMCE thresholds remain fixed. The reference value must also
appear in `BUDGET_VALUES`.

Each group receives `results.csv`. Important fields are `budget`, `wall_ms`,
`budget_fallbacks`, `count_match`, and `slowdown_vs_reference`. Numeric budget
zero forces fallback at the first compatibility examination; `unlimited`
removes the work limit. Optional filters are `BUDGET_DATASETS` and
`BUDGET_BUILD_JOBS`. Reuse the result directory to resume.

### 5.2 Seed-mask capacity ablation

```bash
git switch ablation-capacity
git pull --ff-only origin ablation-capacity

CAPACITY_REPETITIONS=3 \
CAPACITY_TIMEOUT_SECONDS=1800 \
CAPACITY_FAIL_ON_ERROR=1 \
  ./ablation/capacity/run.sh \
    /absolute/path/to/adjacency \
    /absolute/path/to/results/capacity
```

This builds capacities 64, 128, 256, 526, and an unlimited dynamic mask while
holding budget 1,000 and the remaining algorithm configuration fixed. Capacity
128 is run first as the reference. Each group receives raw `results.csv` and a
derived `summary.csv` with the geometric-mean wall-time ratio and speedup versus
128, mean resource use, fallback totals, maximum observed constraint count, and
exact-count agreement. A `geomean_speedup_vs_128` above 1 favors the selected
capacity. Optional filters are `CAPACITY_DATASETS` and `CAPACITY_BUILD_JOBS`.

### 5.3 Retained-clique memory ablation

```bash
git switch ablation-memory
git pull --ff-only origin ablation-memory

MEMORY_INTERVAL_US=1000 \
MEMORY_SYSTEMS=all \
MEMORY_TIMEOUT_SECONDS=1800 \
MEMORY_FAIL_ON_ERROR=1 \
  ./ablation/memory/run.sh \
    /absolute/path/to/adjacency \
    /absolute/path/to/edges \
    /absolute/path/to/results/memory
```

This compares Reorder, HBBMC++, Tomita, and standalone CCRMCE while every
system retains the actual vertex list of every reported clique. The sampler
records `VmRSS` and `VmHWM` from `/proc/PID/status`; conversion to Tomita format
is outside the sampled process. Select a comma-separated subset such as
`MEMORY_SYSTEMS=reorder,tomita` for a shorter run.

Each group receives `results.csv` plus `traces/GRAPH/SYSTEM.txt`. Use
`peak_sampled_rss_kb` as the sampled curve peak, `peak_observed_hwm_kb` as the
largest observed kernel high-water mark, and `wait4_peak_rss_kb` as the final
child-process peak. Check `count_match_reorder` before comparing systems. The
default Reorder budget is 1,000 and can be changed with
`MEMORY_REORDER_BUDGET`.

### 5.4 Recursive-state ablation

```bash
git switch ablation-pxr-states
git pull --ff-only origin ablation-pxr-states

PXR_TIMEOUT_SECONDS=1800 \
PXR_FAIL_ON_ERROR=1 \
  ./ablation/pxr_states/run.sh \
    /absolute/path/to/adjacency \
    /absolute/path/to/edges \
    /absolute/path/to/results/pxr-states
```

The historical branch name is retained even though the final Reorder engine
uses CCRMCE rather than the old PXR implementation. The experiment compares
implementation-native recursive search entries:

- Reorder: FindOne plus exhaustive CCRMCE recursive entries;
- HBBMC: recursive `(R,P,X)` calls with graph reduction and ET disabled;
- Tomita: adjacency-list recursive entries;
- standalone CCRMCE: `BKTomitaRecCCRV3` entries.

These counters exclude parsing, ordering, preprocessing/root setup, Reorder
worklist and hitting-set work, graph reduction, and terminal continuation
expansion. They are implementation-level work measurements, not identical
abstract search-tree nodes. The per-group `results.csv` reports raw state
counts, ratios to Reorder, runtime, invariants, and `all_counts_match`.
Optional filters are `PXR_DATASETS` and `PXR_BUILD_JOBS`.

Return to the production implementation after any ablation:

```bash
git switch main
```

## 6. Reproducibility checklist

For every reported campaign:

1. record `git rev-parse HEAD` and the branch name;
2. keep `environment.txt`, which includes configuration, executable hashes,
   and `uname -a`;
3. retain the raw-source URL, download date, source license, and any ID map;
4. retain the normalized paired graph files or their SHA-256 values;
5. for R-MAT, retain GTgraph version, `N`, requested `M`, seed, and the
   `a,b,c,d` configuration;
6. run competing systems sequentially on the same otherwise idle machine;
7. report timeout and failure rows separately rather than as runtimes;
8. require stored-clique and cross-system clique-count agreement before using
   a measurement.

All experiment build directories, logs, and results are ignored by Git. Use a
new result path for an independent campaign or deliberately reuse a path to
resume it.

## 7. Troubleshooting

- **No valid paired graph inputs matched:** check that both roots contain the
  same relative `GROUP/GRAPH` paths and that each adjacency file starts with
  two nonnegative integers.
- **Input validation fails:** remove self-loops and duplicates, use contiguous
  IDs in `0..n-1`, preserve isolate rows, and ensure the adjacency header counts
  each undirected edge once.
- **A campaign returns success despite a timeout:** the default runners record
  failures and continue. Set the experiment's `*_FAIL_ON_ERROR=1` variable for
  a nonzero final exit status.
- **A rerun skips rows:** result directories are resumable. Choose a new result
  directory for a clean repetition.
- **Illegal-instruction failure:** Release builds use `-march=native`; rebuild
  on the machine that will execute the binaries.
- **Memory ablation has no samples:** confirm Linux `/proc` is mounted and
  readable, and increase `MEMORY_INTERVAL_US` if sampling overhead matters.

## Repository layout

```text
our/                         final Reorder/AdaptiveBK implementation
compare/HBBMCPaperFaithful/  independent HBBMC++ implementation
compare/Tomita/              sparse adjacency-list Tomita implementation
compare/CCRMCE/              standalone CCRMCE implementation
experiments/walltime/        comparison adapters and build definitions
scripts/run_benchmark.sh     main four-system experiment
scripts/prepare_graph.sh     graph normalization helper
data/                        bundled processed graph inputs
```

## License

See [`LICENSE`](LICENSE). Third-party implementations and datasets may carry
their own licenses and citation requirements; consult their source directories
and original download pages.
