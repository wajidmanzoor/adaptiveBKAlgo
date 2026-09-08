# AdaptiveBK / ReorderSib

This repository contains an exact maximal-clique enumerator based on the
optimized ReorderSib worklist algorithm. The implementation uses the final
portable configuration selected by the optimization study; experimental,
ablation, profiling, and PGO build paths are intentionally excluded.

The search combines:

- ascending degeneracy reordering and CSR adjacency storage;
- exact sibling generation through a fixed 128-constraint hitting-set solver;
- a bounded seed solver with exact PXR fallback;
- complete-P, complement-matching, and complement-degree-two terminals;
- compact clique storage and exact duplicate detection; and
- specialized small-branch enumeration for up to four expansion vertices.

## Requirements

- CMake 3.16 or newer
- a C++17 compiler

## Build

```bash
cmake -S our -B our/build -DCMAKE_BUILD_TYPE=Release
cmake --build our/build --parallel
```

The executable is `our/build/adaptive_bk`.

## Input format

Input is a symmetric adjacency list. The first line is `n m`, where `n` is the
number of vertices and `m` is the number of undirected edges. It is followed by
exactly `n` rows in vertex order. Each row begins with its vertex ID and then
lists its neighbors:

```text
4 4
0 1 3
1 0 2
2 1 3
3 0 2
```

Vertex IDs must be in `0..n-1`, self-loops are rejected, every undirected edge
must occur in both adjacency rows, and the header counts each edge once.

## Run

```bash
./our/build/adaptive_bk PATH/TO/GRAPH
```

The default configuration reports maximal cliques of size at least three and
uses a seed-solver work budget of 1,000.

```text
Usage: adaptive_bk GRAPH [options]

  --budget N|unlimited  Compatibility-work budget per seed solver call
  --min-clique-size N   Output threshold (default: 3)
  --print-cliques       Print canonical original vertex IDs
  -h, --help            Show help
```

Use `--min-clique-size 1` for conventional maximal-clique enumeration,
including maximal singleton and two-vertex cliques. The output includes the
exact clique count, stored count, selected threshold and budget, and algorithm
runtime in milliseconds. With `--print-cliques`, each result is emitted as a
sorted line beginning with `clique`.

## Reference comparison

`compare/HBBMCPaperFaithful` contains the independent HBBMC++ reference used
for correctness and runtime comparisons. Given paired HBBMC edge-list and
AdaptiveBK adjacency-list inputs, run:

```bash
bash scripts/run_benchmark.sh DATA_ROOT [RESULT_ROOT]
```

Expected input layout:

```text
DATA_ROOT/
├── hbbmc/  # NAME.edges or NAME.txt.clean
└── pure/   # NAME.graph or NAME.txt
```

Run `bash scripts/run_benchmark.sh --help` for available controls. Results,
commands, logs, clique-count checks, and resource measurements are retained in
the selected result directory.

## Source layout

```text
our/
├── CMakeLists.txt
├── main.cpp
├── inc/
└── src/
```

The production implementation has one supported algorithm configuration. It
does not contain runtime ablation switches, optional RMCE preprocessing,
profiling counters, or compiler profile-generation flags.

## License

See [LICENSE](LICENSE).
