#pragma once

#include "graph.hpp"
#include "truss_order.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace hbbmc_faithful {

// This is the heart of the implementation. High-level shape:
//
//  1. Every isolated vertex (no edges at all) is reported directly as a
//     maximal 1-clique (run()).
//  2. For every other maximal clique (size >= 2), the truss order
//     (truss_order.hpp) is used to pick one canonical "owner" edge: the
//     earliest edge (by truss-peeling rank) contained in the clique. Each
//     edge, visited in truss-rank order, is the root of one independent
//     search (enumerate_root_edge / vertex_bk) that only needs to find
//     maximal cliques it owns -- see the README's "Hybrid branching and
//     ownership" section for the correctness argument (in short: a common
//     neighbor is only ever placed in the candidate set P when the edge
//     connecting it to the branch vertex is later than the owner edge in
//     truss rank, and everything else goes to the excluded set X, which
//     still blocks non-maximal output the normal BK way).
//  3. Inside a root's subproblem, vertex_bk is a standard pivoted
//     Bron-Kerbosch recursion over (partial clique R, candidates P,
//     excluded X), optionally sped up by three orthogonal, independently
//     togglable mechanisms layered on top of the base recursion:
//       - RMCE dynamic reduction (apply_rmce_dynamic_reduction): local P/X
//         simplifications proven safe by the RMCE paper's Lemmas 5/7/8.
//       - RMCE forbidden-set reduction (apply_rmce_forbidden_set_reduction):
//         drops dominated members of X (Lemma 9).
//       - t-plex early termination (try_early_termination): when X is
//         empty and P is a t-plex (t<=3), short-circuits the remaining
//         branch-by-branch recursion with direct enumeration (see
//         early_termination.hpp).
//     All three are only applied when the current candidate graph still
//     matches the *original* graph's edges restricted to owner-later edges
//     (candidate_graph_matches_original) -- once earlier-ranked adjacent
//     vertices have been filtered into X by rank-aware branching, these
//     shortcuts are skipped and plain pivoted BK continues instead, per the
//     README's "rank filter" caveat.
struct EnumerationOptions {
  // 0 gives bare HBBMC. Values 1..3 enable the paper's t-plex ET up to t.
  int early_termination_threshold = 0;
  // These two switches are jointly enabled by the complete "rmce" module.
  bool rmce_dynamic_reduction = false;
  bool rmce_forbidden_set_reduction = false;
  bool validate_invariants = false;
  // Reporting filter only. Search, branching, and ET behavior are unchanged.
  std::size_t minimum_output_clique_size = 1U;
};

struct StructuralCounters {
  std::uint64_t root_edge_branches = 0;
  std::uint64_t descendant_edge_branches = 0;
  std::uint64_t root_candidate_vertices = 0;
  std::uint64_t root_excluded_vertices = 0;
  std::uint64_t root_candidate_max = 0;
  std::uint64_t root_excluded_max = 0;
  std::uint64_t root_support_mismatches = 0;
  std::uint64_t root_tau_bound_violations = 0;
  std::uint64_t root_partition_overlap_violations = 0;

  std::uint64_t vertex_recursive_calls = 0;
  std::uint64_t pivot_calls = 0;
  std::uint64_t pivot_branches = 0;
  std::uint64_t max_vertex_depth = 0;

  std::uint64_t dynamic_reduction_calls = 0;
  std::uint64_t dynamic_degree0_removed = 0;
  std::uint64_t dynamic_degree0_outputs = 0;
  std::uint64_t dynamic_degree1_removed = 0;
  std::uint64_t dynamic_degree1_outputs = 0;
  std::uint64_t dynamic_universal_moved = 0;
  std::uint64_t forbidden_set_removed = 0;

  std::uint64_t et_checks = 0;
  std::uint64_t et1_calls = 0;
  std::uint64_t et2_calls = 0;
  std::uint64_t et3_calls = 0;
  std::uint64_t et1_outputs = 0;
  std::uint64_t et2_outputs = 0;
  std::uint64_t et3_outputs = 0;

  std::uint64_t terminal_outputs = 0;
  std::uint64_t isolated_vertex_outputs = 0;
  std::uint64_t owner_violations = 0;
  std::uint64_t invalid_clique_outputs = 0;
  std::uint64_t nonmaximal_outputs = 0;
  std::uint64_t duplicate_outputs = 0;
};

struct EnumerationResult {
  std::uint64_t maximal_clique_count = 0;
  std::vector<std::vector<int>> cliques;
  TrussOrder truss_order;
  StructuralCounters counters;
};

class Enumerator {
public:
  Enumerator(const Graph &graph, EnumerationOptions options = {});
  EnumerationResult run();

private:
  // Builds the root (R={u,v}, P, X) subproblem for the truss-rank `rank`
  // edge `edge_id` and hands it to vertex_bk. This is where the
  // candidate/excluded split by "is the triangle-closing edge later than
  // the owner" (the ownership invariant) is first established.
  void enumerate_root_edge(int edge_id, int rank);
  // The pivoted BK recursion itself: one call handles one (R, P, X) node.
  // `owner_edge` is threaded through unchanged for the whole subtree so
  // every recursive call can keep checking "is this still later than my
  // root's owner edge".
  void vertex_bk(std::vector<int> &partial, std::vector<int> candidates,
                 std::vector<int> excluded, int owner_edge,
                 std::uint64_t depth);
  // RMCE Lemmas 5/7/8 (dynamic degree-0/1/universal reductions). Mutates
  // `candidates`/`excluded` in place and may push forced vertices onto
  // `partial`. Returns true if the reduction proved the current `partial`
  // can no longer lead to a valid (maximal, non-duplicate) output and the
  // caller should stop without emitting it.
  bool apply_rmce_dynamic_reduction(std::vector<int> &partial,
                                    std::vector<int> &candidates,
                                    std::vector<int> &excluded, int owner_edge);
  // RMCE Lemma 9 (forbidden-set containment): drops members of `excluded`
  // whose forbidden neighborhood is dominated by another member's.
  void apply_rmce_forbidden_set_reduction(const std::vector<int> &candidates,
                                          std::vector<int> &excluded);
  // If X is empty and P is a t-plex for some t <= options_.early_termination_threshold,
  // enumerates every maximal clique of G[P] directly (early_termination.hpp)
  // instead of continuing to branch, and returns true. Returns false (no-op)
  // otherwise.
  bool try_early_termination(const std::vector<int> &partial,
                             const std::vector<int> &candidates,
                             const std::vector<int> &excluded, int owner_edge);
  // Picks the candidate-or-excluded vertex with the most neighbors inside P
  // (classic BK pivoting: branching only on P \ N(pivot) is what keeps the
  // recursion from redundantly exploring cliques that share a big common
  // neighborhood).
  int choose_pivot(const std::vector<int> &candidates,
                   const std::vector<int> &excluded) const;
  // True iff edge (u,v) has a strictly later truss rank than `owner_edge`
  // (or, when owner_edge==-1, i.e. no root context, just falls back to plain
  // adjacency). This is the O(1) primitive the whole ownership scheme is
  // built on.
  bool edge_is_after_owner(int u, int v, int owner_edge) const;
  // True iff every edge among `candidates` that exists in the original
  // graph is also later than `owner_edge` -- i.e. the P/X split at this
  // node still exactly mirrors the original induced subgraph, so the RMCE
  // and ET shortcuts (which assume they're operating on an unmodified
  // induced subgraph) remain sound to apply.
  bool candidate_graph_matches_original(const std::vector<int> &candidates,
                                        int owner_edge) const;
  // The single output sink used by every terminal, ET, and reduction path.
  // It always increments the count and retains the complete clique identity.
  void emit_clique(std::vector<int> clique, int owner_edge);

  bool is_clique(const std::vector<int> &clique) const;
  bool is_maximal(const std::vector<int> &clique) const;
  std::string clique_key(const std::vector<int> &clique) const;

  const Graph &graph_;
  EnumerationOptions options_;
  EnumerationResult result_;
  // Validation-mode-only duplicate detector, keyed by clique_key().
  std::unordered_set<std::string> seen_cliques_;
  // Epoch-stamped "is this vertex adjacent to something in X" marks used by
  // apply_rmce_dynamic_reduction's relaxed degree-one rule (Lemma 7); reset
  // in O(1) per call by bumping dynamic_mark_epoch_value_ instead of
  // clearing the whole vector (see that function for the wraparound guard).
  std::vector<std::uint32_t> dynamic_mark_epoch_;
  std::uint32_t dynamic_mark_epoch_value_ = 0U;
};

} // namespace hbbmc_faithful
