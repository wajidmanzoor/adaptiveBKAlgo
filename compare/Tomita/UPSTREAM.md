# Tomita maximal-clique implementation: upstream snapshot

This directory started from the tracked files in **Quick Cliques v1.0**,
downloaded from:

- Repository: https://github.com/darrenstrash/quick-cliques
- Tag: `v1.0`
- Commit: `dd9c4d9fdf8243e2861aa6b99e8992836da741a9`
- Retrieved: 2026-09-07

The Tomita implementation is in:

- `src/tomita.c` (program entry point)
- `src/tomita_algorithm.c` (algorithm implementation)
- `src/tomita_algorithm.h` (public declarations)

This is Darren Strash's public C implementation of the maximal-clique
enumeration algorithm described by Tomita, Tanaka, and Takahashi (TCS 2006).
It should not be described as source code published by the paper's authors.
The upstream license is GNU GPL v3; see `COPYING` and `LICENSE`.

## Local benchmark changes

The local copy is configured to retain every maximal clique of size at least
three in memory and not print clique identities. The algorithm timer uses
`CLOCK_MONOTONIC` and includes insertion of each retained identity. Aggregate
configuration, count, storage, and timing keys are still printed. The four
upstream enumeration implementations use the same minimum-size and storage
policy, although the twelve-dataset runner uses only the sparse
`tomita-adjacency-list` implementation.

The original `tomita` executable uses an `n` by `n` byte adjacency matrix and
is therefore not viable for the largest sparse datasets. Quick Cliques also
ships `tomita-adjacency-list`, which preserves Tomita's pivoted search while
using linear graph storage. Results from it must be labeled
`tomita-adjacency-list`, not matrix Tomita.

## Build on current GCC

The original Makefile predates modern C inline semantics. The local Makefile
now enables the intended GNU89 inline behavior, so build with:

```sh
make
```

Run the Tomita implementation with:

```sh
./bin/tomita GRAPH
```

Reading from standard input (`./bin/tomita < GRAPH`) remains supported. The
same two forms work with `./bin/adjlist`.

The original input format has the vertex count on line 1, twice the undirected
edge count on line 2, and then both directed forms of every edge as `u,v`.
