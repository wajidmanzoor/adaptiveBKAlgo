# Independent HBBMC/HBBMC++ correctness reference

This directory is a clean-room, from-paper implementation of maximal-clique
enumeration (MCE) as described in *Maximal Clique Enumeration with Hybrid
Branching and Early Termination* (arXiv:2412.08218) and of the graph-reduction
rules from *Accelerating Maximal Clique Enumeration via Graph Reduction*
(PVLDB 17(10), 2024). It is intended to recover the algorithm that the HBBMC
paper evaluates, with an auditable correctness path and an honest provenance
record.

It is **algorithm- and correctness-faithful, not performance-equivalent to the
authors' unavailable experimental executable**. In particular, the graph
representation, truss peeling, and RMCE work queues use straightforward STL
containers. They preserve the paper's search space and reductions but do not
reproduce undocumented data layouts, low-level optimizations, or compiler
choices. Timings from this binary must be described as timings of this
independent implementation, never as timings of the authors' artifact.

## Why an independent implementation is necessary

The public HBBMC repository was inspected at commit
`baf7f57e8d30c499b394fe3e1cf550dabfb73ac0`. Its CMake target compiles the
utility and `enum.cpp` path; that main invokes `Tomita_opt_mat`, ignores the
second positional option, and the checked-in source does not build unchanged.
It therefore does not expose the paper's HBBMC++ MCE executable.

Repository history was also audited. The historical `HBBMC.zip` present at
commit `370309a` and deleted at `8075894` contains the same edge-oriented source
that is omitted from CMake. Its `EBBMC_plus_plus` path uses a global fixed `K`
and implements fixed-k clique listing, not the paper's hybrid MCE algorithm.
Thus neither current nor historical public source provides a runnable
performance-matched HBBMC++ baseline.

The RMCE reference checkout used to trace its reductions is pinned at commit
`9eca92678a69838fbba29e9279be2535d9929d1e`. No code was copied from either
upstream; the implementation here follows the published invariants and was
differentially validated.

## Implemented variants

| CLI configuration | Reported variant | Active mechanisms |
|---|---|---|
| `--graph-reduction none --et 0` | HBBMC | truss-ordered edge roots; pivoted vertex-BK descendants |
| `--graph-reduction none --et t`, `t=1..3` | HBBMC+ET(t) | HBBMC plus t-plex early termination (diagnostic ablation) |
| `--graph-reduction rmce --et 0` | HBBMC+ | HBBMC plus complete RMCE global, dynamic, and forbidden-set reduction |
| `--graph-reduction rmce --et 1` or `2` | HBBMC+GR+ET(t) | diagnostic ablation with all reductions and the selected ET threshold |
| `--graph-reduction rmce --et 3` | HBBMC++ | the paper's full combination: HBBMC + GR + ET with `t=3` |

The executable prints `hbbmc_plus_plus=true` only for the last configuration,
when global, dynamic, and forbidden-set reductions and ET(3) are all active.
`--early-termination` is an alias for `--et`.

## Feature and proof mapping

### Hybrid branching and ownership

The implementation greedily removes a current minimum-support edge to form the
truss edge order. Every maximal clique of size at least two is assigned to its
earliest edge in that total order. For root edge `e=(u,v)`, a common neighbor
is placed in `P` exactly when both incident triangle edges occur after `e`; all
other common neighbors are placed in `X`. At a descendant branch on `v`, an
adjacent candidate `w` remains in `P` only when edge `(v,w)` also occurs after
the owner edge; an earlier-ranked adjacent vertex is promoted to `X`, where it
continues to block non-maximal output in the original graph. Inductively, every
non-root internal edge of a reported clique is later than its owner, while `X`
retains every original-graph maximality blocker. Consequently, the clique is
present in exactly the root owned by its earliest edge, including when an
earlier internal edge is disjoint from the current branch edge. Descendants use
maximum-candidate-degree pivoted vertex BK. Original isolates are emitted
separately as maximal one-vertex cliques.

### Early termination

For threshold `t<=3`, ET is used only when `X` is empty and `G[P]` is a t-plex.
For `t=1`, `P` is one clique. For `t=2`, the complement of `G[P]` is a matching.
For `t=3`, every complement component is an isolated vertex, path, or cycle;
the path/cycle recurrences enumerate exactly its maximal independent sets,
which are exactly the maximal cliques of `G[P]`. The enumerator visits every
terminal continuation and constructs `R` plus its clique suffix, including in
count-only mode. Count-only execution discards that identity after incrementing
the output count; it does not replace the traversal with an algebraic aggregate.

### RMCE reductions

The final-paper numbering mapped to the code is:

| RMCE result | Implementation | Preservation argument |
|---|---|---|
| Lemmas 1--3, global vertex reduction (Algorithm 3) | iterative degree-0/1/2 queue in `reduction.cpp` | a degree-1 edge is maximal; a degree-2 vertex yields either two maximal edges or one maximal triangle; all incident search space is then deleted |
| Lemma 4, global edge reduction (Algorithm 4) | repeated deletion of edges in no triangle | such an edge cannot be extended, so it is itself a maximal clique |
| Lemma 5, dynamic degree zero | candidate-local degree-zero removal in `enumerator.cpp` | emit `R+u` iff no vertex in `X` extends it; otherwise remove it silently |
| Lemma 7, relaxed dynamic degree one | marked-`P` scan using adjacency to `X` | if either endpoint has no neighbor in `X`, no forbidden vertex extends both, so `R+u+v` is maximal |
| Lemma 8, dynamic degree `|P|-1` | move universal candidates into `R` and intersect `X` | every maximal clique in the subproblem must contain a candidate adjacent to all other candidates |
| Lemma 9, forbidden-set containment (Algorithm 6) | remove `x_i` when `N_P(x_i)` is contained in another forbidden neighborhood | every suffix blocked by the removed vertex is also blocked by the dominating forbidden vertex |

The global degree-zero rule is adapted to the HBBMC paper's definition of MCE:
an isolate present in the input is emitted as a maximal singleton. A vertex
that becomes isolated because prior rules consumed all its incident search
space is removed silently. In the adjacent degree-two case, the surviving
neighbor edge is deleted only when the reduced vertex is its sole common
neighbor; otherwise it remains available for a different maximal clique.
Global vertex and edge phases repeat to a fixed point because edge deletion can
create new low-degree vertices.

Every global output contains a structure that is deleted immediately, while
the residual graph retains only unconsumed search space. Dynamic reductions
apply the corresponding `P/X` identities inside each edge-owned BK subproblem.
Together these facts make direct and residual outputs disjoint and preserve the
complete MCE set. Validation mode checks this property explicitly.

The rank filter can make the owner-specific candidate graph a strict subgraph
of the original induced graph. RMCE dynamic/forbidden reductions and plex early
termination are therefore invoked only when every original edge among the
current candidates is later than the owner edge; otherwise exact rank-aware BK
recursion continues. Each recursive call also restores the mutable partial
clique to its entry size, so a universal candidate moved into `R` by RMCE is
local to that call and cannot leak into sibling branches.

The forbidden-set implementation evaluates Lemma 9 containment directly inside
each edge-owned root subproblem. It does not reproduce RMCE's degeneracy-root
`ignoreId` cache, whose layout is tied to a vertex-root search. This preserves
the lemma's semantics but is another deliberate engineering/performance
difference from the author artifacts.

## Build and run

```sh
cmake -S compare/HBBMC-faithful -B /tmp/hbbmc-faithful-build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/hbbmc-faithful-build -j

/tmp/hbbmc-faithful-build/hbbmc_faithful GRAPH \
  --graph-reduction rmce --early-termination 3
```

Input is a whitespace-separated undirected edge list. Blank lines and lines
starting with `#` or `%` are ignored; loops and duplicate edges are removed.
Use `--num-vertices N` when labels are `0..N-1` and isolated vertices must be
retained. Without it, only labels appearing in an edge are observable.

By default the executable is count-only. It does not store output identities,
does not allocate duplicate-detection sets, skips clique/maximality/ownership
validation, but it still enumerates and temporarily constructs every clique
reached through early termination. `--print-cliques` retains those identities.
`--validate` additionally checks clique validity, global maximality, uniqueness,
and earliest-edge ownership. Reported `algorithm_runtime_ms` includes graph
reduction, truss ordering, and enumeration, but excludes input parsing; the
component times are also printed.

## Validation evidence

`tests/test_exhaustive.cpp` performs the following checks:

- exhaustive reference comparison on all 33,868 labeled simple graphs with
  `0<=n<=6`;
- 135,472 unreduced runs covering ET thresholds 0, 1, 2, and 3;
- 67,736 complete RMCE pipelines, with both HBBMC+ (`t=0`) and HBBMC++ (`t=3`)
  compared against the unreduced maximal-clique set;
- count-only versus materialized/validated count equality on every reduced
  pipeline, while count-only result identity vectors remain empty;
- direct ET path and cycle fixtures through eight vertices and an integrated
  15-vertex HBBMC++ fixture that reaches the actual ET3 sink after all RMCE
  reduction paths are enabled;
- a 16-vertex counterexample with a disjoint earlier internal edge, which
  verifies unique earliest-edge ownership, and a 262-vertex fixture that forces
  RMCE universal moves before sibling branches, which verifies recursive
  partial-state restoration;
- structural checks for root ownership, truss support bounds, duplicates,
  invalid/non-maximal outputs, and coverage of every global and dynamic rule.

Run the suite with:

```sh
ctest --test-dir /tmp/hbbmc-faithful-build --output-on-failure
```

Fresh Release and ASan/UBSan builds pass the full suite. Corrected count-only
HBBMC++ runs also reconcile the public totals on DBLP (125,568 direct plus
131,983 residual equals 257,551) and YouTube (1,688,133 direct plus 1,577,823
residual equals 3,265,956). These larger-graph checks are count reconciliation,
not exhaustive materialized-identity proofs. The bare-HBBMC lane independently
reports the same two totals, isolating the ownership repair from RMCE and plex
termination. Together the evidence establishes a strong correctness result;
it is not a proof of performance equivalence.

## Performance assessment and fair use

The count-only hot path avoids persistent output-identity and duplicate-set
storage, and owner checks remain validation-only. RMCE direct outputs are still
counted without identity storage, while ET regions explicitly enumerate every
suffix and temporarily construct every resulting clique. This remains a
**correctness reference, not a performance-grade reproduction of the authors'
system**. The set-based global reducer can have a substantially larger constant
factor (and a weaker engineering bound) than RMCE's specialized arrays, and the
truss implementation uses hash sets and a lazy heap. Any paper using its runtime
should label it “independent paper-faithful HBBMC++ reimplementation,” disclose
this limitation, compile all methods comparably, include preprocessing for every
method, and avoid attributing observed runtime differences to the original
HBBMC authors' code.
