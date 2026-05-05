---
name: Project Overview
description: Adaptive BK maximal clique enumeration research — class hierarchy, key algorithms, and design intent
type: project
---

C++ research project implementing variants of Bron-Kerbosch (BK) maximal clique enumeration. All code lives in `orginal/` (the "original" branch of experiments). Main classes in `orginal/src/helpers.cpp` / `orginal/inc/helpers.h`.

**Why:** Research into adaptive pruning strategies that reuse information from already-found cliques to prune future branches.

**How to apply:** Suggestions should preserve the `mustin`/`expandTo` framework and the sibling-set hitting-set design.

## Class Hierarchy (modes in main.cpp)

- `AdjMatBK` (modes unused directly) — baseline BK on adjacency matrix
- `AdjListBK` (modes unused directly) — baseline BK on adjacency list
- `PivotBK` (modes 0–2) — standard BK with max-neighborhood pivot; DegOrder controls vertex ordering
- `Reorder` (modes 3–5) — custom mustin/expandTo framework; after finding clique C, re-seeds sibling branches
- `ReorderSib` (modes 6–10) — extends Reorder with explicit sibling-set computation via hitting-set solvers
- `ReorderSibOwner` (mode 11) — declared but NOT YET IMPLEMENTED in helpers.cpp

## Core Data Structures (ReorderSib)

- `mustin[i]`: mandatory vertex set for branch i (must appear in the clique)
- `expandTo[i]`: candidate vertices for branch i
- `allCliques`, `foundLevel`, `cliquesByVertex`: track discovered cliques per vertex and level
- `fullSkipCheck[i]`: if true, branch i's entire search space (M ∪ E) is checked against new cliques for aggressive pruning
- `seenCliques`: deduplication set (string-encoded clique keys)

## Key Algorithms

**collectCoveringCliques(M, level):** Finds already-found cliques containing M. Uses the vertex with fewest cliques as seed for efficiency.

**buildHitSets(E, cliqueIds):** For each covering clique C, constraint = E \ C (must pick something outside C).

**generateSiblingSetsFromCliques(E, cliqueIds):** Solves hitting-set problem over hit sets. Strategy selected by `SibMethod` enum: BRUTE_FORCE, BACKTRACKING, GREEDY, BITMASK, MIN_HITTING_SET.

**minimalByInclusion(solutions):** Removes dominated sets — if S1 ⊇ S2, drop S1.

**Reorder update (after finding clique C):** `expandTo[i] ← intersect(setDiff(adjList[mustin[i].back()], mustin[i]), unionSet(C, expandTo[i]))`. Allows branches to "see" C as a new expansion region.

**fullSkipCheck pruning:** When a branch was seeded with covering-clique info, skip it entirely if `M ∪ E ⊆ C` for the newly found clique C.

## Build

`orginal/` has a `Makefile` and `CMakeLists.txt`. Binary: `orginal/bk_algorithm`. Usage: `./bk_algorithm <graph_file> <mode>` where mode is 0–11.
