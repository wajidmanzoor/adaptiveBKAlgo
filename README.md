# PivotBK and Pure ReorderSib

This repository provides two maximal-clique enumeration modes:

- mode `0`: PivotBK, the conventional pivoting baseline;
- mode `1`: Pure ReorderSib, the optimized exact worklist implementation used
  by the comparison benchmark.

It also includes an independent paper-faithful HBBMC++ implementation used as
the reference.

## Requirements

- A C++17 compiler
- CMake 3.16 or newer
- Linux/Ubuntu with GNU `time`, `timeout`, `sha256sum`, and `strip` for the
  comparison benchmark

All commands below are run from the repository root.

## Build

```bash
cmake -S our -B our/build -DCMAKE_BUILD_TYPE=Release
cmake --build our/build --parallel
```

This creates `our/build/bk_algorithm`.

## Input format

The main executable expects an adjacency list. The first line is `n m`; each
following line starts with a vertex ID and lists its neighbors:

```text
4 4
0 1 3
1 0 2
2 1 3
3 0 2
```

Vertex IDs must be in `0..n-1`. Each undirected edge must appear in both
vertices' adjacency lists, while `m` counts each edge once.

## Run mode 0: PivotBK

```bash
./our/build/bk_algorithm PATH/TO/GRAPH 0 0
```

The final method argument is required by the common CLI but is ignored by
PivotBK. Both algorithms always use ascending degeneracy order. The optional
next argument sets the minimum reported clique size:

```bash
# Include singleton and two-vertex maximal cliques.
./our/build/bk_algorithm PATH/TO/GRAPH 0 0 1

# Report only maximal cliques of size at least four.
./our/build/bk_algorithm PATH/TO/GRAPH 0 0 4
```

The default minimum clique size is `3`.

## Run mode 1: Pure ReorderSib

The complete command is:

```bash
./our/build/bk_algorithm GRAPH MODE METHOD \
  [hitSetLimit] \
  [prune1] [prune2] \
  [sp1] [sp2] [sp3] [sp4] [sp5] [sp6] \
  [minCliqueSize]
```

The standard optimized configuration is:

```bash
cmake -S our -B our/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DREORDERSIB_PROFILING=OFF \
  -DPURE_HITSET_VARIANT=128
cmake --build our/build --parallel

PURE_HITSET_BUDGET=10000 \
  ./our/build/bk_algorithm PATH/TO/GRAPH \
  1 1 4294967295 1 1 1 1 1 1 1 1 3
```

This selects optimized exact hitting-set search with the fixed ascending
degeneracy order, all solver flags, minimum clique size three, a 128-entry
hitting-set build, and a 10,000-work solver budget. The Pure exact lane retains
the legacy pruning-flag positions for command compatibility; `sp1`, `sp2`, and
`sp3` control its hitting-set solver.

Use a final value of `1` to include singleton and two-vertex maximal cliques.
Set `PURE_RMCE=1` to enable optional RMCE preprocessing.

## Run the paper-faithful HBBMC++ reference directly

The reference accepts a whitespace-separated undirected edge list with one
edge `u v` per line:

```bash
cmake -S compare/HBBMCPaperFaithful \
  -B compare/HBBMCPaperFaithful/build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build compare/HBBMCPaperFaithful/build --parallel

./compare/HBBMCPaperFaithful/build/hbbmc_faithful PATH/TO/EDGES \
  --graph-reduction rmce \
  --et 3 \
  --num-vertices N \
  --min-clique-size 3
```

Replace `N` with the graph's vertex count. `--num-vertices N` preserves
isolated vertices whose labels do not occur in the edge list.

## Run Pure ReorderSib vs. the paper-faithful reference

The benchmark requires paired versions of each dataset:

```text
DATA_ROOT/
├── hbbmc/   # edge-list files: NAME.edges or NAME.txt.clean
└── pure/   # adjacency-list files: NAME.graph or NAME.txt
```

Run the comparison on Linux/Ubuntu:

```bash
bash scripts/run_benchmark.sh DATA_ROOT [RESULT_ROOT]
```

For example:

```bash
COMPARE_DATASETS="dblp,youtube" \
COMPARE_TIMEOUT_SECONDS=1200 \
  bash scripts/run_benchmark.sh \
  /data/labdata/wajid/hbbmcData results/reference_run
```

If `RESULT_ROOT` is omitted, the script creates a timestamped directory under
`results/`. It contains `runs.csv`, `comparison.csv`, and `environment.txt`.
Build logs, commands, stdout/stderr, and resource measurements are written to
`logs/` unless `COMPARE_LOG_DIR` overrides that location. Reusing an explicit
result directory resumes a partial campaign and skips recorded rows.

Show every benchmark option without starting a run:

```bash
bash scripts/run_benchmark.sh --help
```
