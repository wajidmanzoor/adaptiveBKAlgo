#include "../inc/enumerator.hpp"

#include "../inc/early_termination.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace hbbmc_faithful {

// All of `candidates`, `excluded`, and neighbor lists (Graph::neighbors) are
// maintained as sorted, duplicate-free vectors throughout this file, which
// is what lets every set operation below run as a linear merge/binary
// search instead of falling back to hash sets.
namespace {

// vertices ∩ neighbors, both assumed sorted.
std::vector<int> intersect_neighbors(const std::vector<int> &vertices,
                                     const std::vector<int> &neighbors) {
  std::vector<int> output;
  output.reserve(std::min(vertices.size(), neighbors.size()));
  std::set_intersection(vertices.begin(), vertices.end(), neighbors.begin(),
                        neighbors.end(), std::back_inserter(output));
  return output;
}

// vertices \ neighbors, both assumed sorted. Used to compute the pivot's
// non-neighbors within P: the classic pivoted-BK branch set.
std::vector<int> subtract_neighbors(const std::vector<int> &vertices,
                                    const std::vector<int> &neighbors) {
  std::vector<int> output;
  output.reserve(vertices.size());
  std::set_difference(vertices.begin(), vertices.end(), neighbors.begin(),
                      neighbors.end(), std::back_inserter(output));
  return output;
}

// |a ∩ b| without materializing the intersection, both assumed sorted.
std::size_t intersection_size(const std::vector<int> &a,
                              const std::vector<int> &b) {
  std::size_t i = 0;
  std::size_t j = 0;
  std::size_t count = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] < b[j]) {
      ++i;
    } else if (b[j] < a[i]) {
      ++j;
    } else {
      ++count;
      ++i;
      ++j;
    }
  }
  return count;
}

// Inserts `value` into the sorted vector `values` if not already present.
void insert_sorted_unique(std::vector<int> &values, int value) {
  const auto position = std::lower_bound(values.begin(), values.end(), value);
  if (position == values.end() || *position != value) {
    values.insert(position, value);
  }
}

// Removes `value` from the sorted vector `values` if present; returns
// whether it was found and removed.
bool erase_sorted(std::vector<int> &values, int value) {
  const auto position = std::lower_bound(values.begin(), values.end(), value);
  if (position == values.end() || *position != value) {
    return false;
  }
  values.erase(position);
  return true;
}

} // namespace

Enumerator::Enumerator(const Graph &graph, EnumerationOptions options)
    : graph_(graph), options_(options) {
  if (options_.early_termination_threshold < 0 ||
      options_.early_termination_threshold > 3) {
    throw std::invalid_argument("early-termination threshold must be in [0,3]");
  }
  if (options_.minimum_output_clique_size == 0U) {
    throw std::invalid_argument(
        "minimum output clique size must be at least one");
  }
}

EnumerationResult Enumerator::run() {
  result_ = {};
  seen_cliques_.clear();
  dynamic_mark_epoch_.assign(static_cast<std::size_t>(graph_.vertex_count()),
                             0U);
  dynamic_mark_epoch_value_ = 0U;
  // Compute the total truss-peeling edge order once up front; every root
  // subproblem below is defined relative to this fixed order (see
  // enumerate_root_edge / edge_is_after_owner).
  result_.truss_order = compute_truss_order(graph_);

  // Vertices with no edges at all can't be an owner-edge's clique member,
  // so they're reported separately here as their own maximal 1-cliques
  // (owner_edge=-1 signals "no owning edge" to emit_clique's validation).
  for (int v = 0; v < graph_.vertex_count(); ++v) {
    if (graph_.neighbors(v).empty()) {
      ++result_.counters.isolated_vertex_outputs;
      if (options_.collect_cliques || options_.validate_invariants) {
        emit_clique({v}, -1);
      } else {
        record_output_count(1U);
      }
    }
  }

  // Visit every edge, root-first, in ascending truss rank. Rank order
  // matters here only insofar as each root needs to know its own rank (to
  // compare triangle edges against); roots are otherwise independent.
  for (std::size_t rank = 0; rank < result_.truss_order.edge_at_rank.size();
       ++rank) {
    enumerate_root_edge(result_.truss_order.edge_at_rank[rank],
                        static_cast<int>(rank));
  }
  return std::move(result_);
}

// Sets up the (R={u,v}, P, X) root subproblem owned by edge (u,v) at truss
// rank `rank`, then delegates to vertex_bk. For every common neighbor w of
// u and v (i.e. every vertex that closes a triangle u-v-w), w goes into P
// only if *both* triangle-closing edges (u,w) and (v,w) rank after this
// root edge; otherwise w goes into X. This is the ownership invariant: any
// maximal clique built from here can only grow through vertices connected
// via later edges, which is exactly what guarantees the clique is later
// found to be *owned* by this edge and not some earlier one (see the
// class-level comment in enumerator.hpp and the README's proof sketch).
void Enumerator::enumerate_root_edge(int edge_id, int rank) {
  ++result_.counters.root_edge_branches;
  const Edge edge = graph_.edge(edge_id);
  std::vector<int> candidates;
  std::vector<int> excluded;
  for (const int w : graph_.common_neighbors(edge.u, edge.v)) {
    const int uw = graph_.edge_id(edge.u, w);
    const int vw = graph_.edge_id(edge.v, w);
    if (uw < 0 || vw < 0) {
      throw std::logic_error("triangle edge is absent from edge map");
    }
    if (result_.truss_order.rank_of_edge[static_cast<std::size_t>(uw)] > rank &&
        result_.truss_order.rank_of_edge[static_cast<std::size_t>(vw)] > rank) {
      candidates.push_back(w);
    } else {
      excluded.push_back(w);
    }
  }
  std::sort(candidates.begin(), candidates.end());
  std::sort(excluded.begin(), excluded.end());

  result_.counters.root_candidate_vertices += candidates.size();
  result_.counters.root_excluded_vertices += excluded.size();
  result_.counters.root_candidate_max = std::max<std::uint64_t>(
      result_.counters.root_candidate_max, candidates.size());
  result_.counters.root_excluded_max = std::max<std::uint64_t>(
      result_.counters.root_excluded_max, excluded.size());

  // Diagnostics only (not used for control flow): the edge's support at
  // removal time (number of w with both uw,vw not yet peeled, i.e. ranked
  // after this edge) is exactly the same condition used to build
  // `candidates` above, so |P| should always equal support_at_removal
  // exactly, and by definition of tau (the max support seen during
  // peeling) should never exceed it either.
  if (candidates.size() !=
      static_cast<std::size_t>(
          result_.truss_order
              .support_at_removal[static_cast<std::size_t>(edge_id)])) {
    ++result_.counters.root_support_mismatches;
  }
  if (candidates.size() > static_cast<std::size_t>(result_.truss_order.tau)) {
    ++result_.counters.root_tau_bound_violations;
  }
  if (options_.validate_invariants) {
    std::vector<int> overlap;
    std::set_intersection(candidates.begin(), candidates.end(),
                          excluded.begin(), excluded.end(),
                          std::back_inserter(overlap));
    if (!overlap.empty()) {
      ++result_.counters.root_partition_overlap_violations;
    }
  }

  std::vector<int> partial{edge.u, edge.v};
  vertex_bk(partial, std::move(candidates), std::move(excluded), edge_id, 0U);
}

// One node of the pivoted BK recursion tree. `partial` is R, passed by
// reference and mutated in place (pushed/popped around recursive calls)
// rather than copied, since it's the same across all siblings; `candidates`
// (P) and `excluded` (X) are passed by value because each branch needs its
// own independent copy to mutate.
void Enumerator::vertex_bk(std::vector<int> &partial,
                           std::vector<int> candidates,
                           std::vector<int> excluded, int owner_edge,
                           std::uint64_t depth) {
  // Because `partial` is shared/mutated across the recursion instead of
  // copied, RMCE dynamic reduction below can push extra "universal
  // candidate" vertices onto it (see apply_rmce_dynamic_reduction's Lemma-8
  // handling). This guarantees partial is always restored to exactly its
  // entry size on every return path (including exceptions), so a universal
  // vertex moved into R by one branch never leaks into a sibling branch.
  const std::size_t partial_size_on_entry = partial.size();
  struct PartialRestore {
    std::vector<int> &partial;
    std::size_t size;
    ~PartialRestore() { partial.resize(size); }
  } partial_restore{partial, partial_size_on_entry};

  ++result_.counters.vertex_recursive_calls;
  result_.counters.max_vertex_depth =
      std::max(result_.counters.max_vertex_depth, depth);

  // The RMCE reductions and t-plex ET all assume they're looking at a
  // genuine induced subgraph (P, X exactly as they'd be in an unfiltered
  // BK run), which only holds while rank-aware branching hasn't yet
  // demoted an original-graph edge's endpoint from P to X (see
  // candidate_graph_matches_original's doc comment). Once that has
  // happened at some ancestor call, these three mechanisms are skipped for
  // the rest of this subtree and plain pivoted BK takes over.
  //
  // RMCE applies maximality-check reduction before its recursive kernel and
  // dynamic reduction at kernel entry. In HBBMC, each truss-owned edge child
  // is the analogous top-level subproblem.
  const bool owner_graph_unfiltered =
      candidate_graph_matches_original(candidates, owner_edge);
  // Forbidden-set reduction (Lemma 9) only runs once, at the root of each
  // owner's subtree (depth==0), matching RMCE's degeneracy-root scoping
  // (see the README's note on why this isn't a per-node cache).
  if (options_.rmce_forbidden_set_reduction && depth == 0U &&
      owner_graph_unfiltered) {
    apply_rmce_forbidden_set_reduction(candidates, excluded);
  }
  bool partial_blocked_by_reduced_candidate = false;
  if (options_.rmce_dynamic_reduction && owner_graph_unfiltered) {
    partial_blocked_by_reduced_candidate =
        apply_rmce_dynamic_reduction(partial, candidates, excluded, owner_edge);
  }

  // Standard BK base case: P and X both empty means `partial` (R) is a
  // maximal clique (X empty proves no excluded vertex could extend it,
  // which is what "maximal" requires). If dynamic reduction discarded a
  // candidate that could still have extended partial without proving it
  // unreachable another way, `partial_blocked_by_reduced_candidate`
  // suppresses this output since it would otherwise be non-maximal.
  if (candidates.empty()) {
    if (excluded.empty() && !partial_blocked_by_reduced_candidate) {
      ++result_.counters.terminal_outputs;
      if (options_.collect_cliques || options_.validate_invariants) {
        emit_clique(partial, owner_edge);
      } else {
        record_output_count(partial.size());
      }
    }
    return;
  }

  // Early termination short-circuit: only tried when P is nonempty (handled
  // above) and the induced-subgraph precondition still holds.
  if (owner_graph_unfiltered &&
      try_early_termination(partial, candidates, excluded, owner_edge)) {
    return;
  }

  // Pivoted BK branching: only vertices in P that are NOT neighbors of the
  // pivot need their own branch, since any maximal clique extension through
  // a P-neighbor of the pivot is already covered by the branch that
  // eventually selects the pivot itself (or another P-neighbor of it).
  const int pivot = choose_pivot(candidates, excluded);
  if (pivot < 0) {
    throw std::logic_error("nonempty candidate set has no pivot");
  }
  ++result_.counters.pivot_calls;
  const std::vector<int> branch_vertices =
      subtract_neighbors(candidates, graph_.neighbors(pivot));

  // branch_vertices is a fixed snapshot taken before the loop starts, while
  // `candidates` itself shrinks by one (the just-visited v) at the bottom
  // of every iteration; since branch_vertices has no duplicate entries,
  // every subsequent v is still present in `candidates` at the top of its
  // own iteration, so this lookup is a defensive check rather than a
  // reachable skip in current usage.
  for (const int v : branch_vertices) {
    const auto position =
        std::lower_bound(candidates.begin(), candidates.end(), v);
    if (position == candidates.end() || *position != v) {
      continue;
    }
    ++result_.counters.pivot_branches;
    // Standard BK move-to-R step: descend with v added to R, and build the
    // child's P/X by intersecting with v's neighborhood -- both X (via
    // next_excluded, straightforward set intersection) and P (via the loop
    // below, which additionally re-applies the ownership rank check).
    std::vector<int> next_candidates;
    next_candidates.reserve(candidates.size());
    std::vector<int> next_excluded =
        intersect_neighbors(excluded, graph_.neighbors(v));
    // A root owns only cliques whose remaining internal edges occur after
    // it in the truss order. An adjacent but earlier-ranked vertex remains
    // an original-graph maximality blocker, so promote it from P to X.
    for (const int w : candidates) {
      if (w == v || !graph_.adjacent(v, w)) {
        continue;
      }
      if (edge_is_after_owner(v, w, owner_edge)) {
        next_candidates.push_back(w);
      } else {
        insert_sorted_unique(next_excluded, w);
      }
    }
    partial.push_back(v);
    vertex_bk(partial, std::move(next_candidates), std::move(next_excluded),
              owner_edge, depth + 1U);
    partial.pop_back();

    // Standard BK move-to-X step for the *sibling* branches still to come:
    // v is done contributing new branches, so remove it from P and add it
    // to X (any clique through v was already fully explored above).
    candidates.erase(std::lower_bound(candidates.begin(), candidates.end(), v));
    insert_sorted_unique(excluded, v);
  }
}

// Applies three RMCE "dynamic" reductions (so called because they operate
// on the current P/X of one BK node, unlike the global vertex/edge rules in
// reduction.cpp which run once on the whole graph up front) in sequence:
//   1. Lemma 5, degree-0 in P: a candidate v with no P-neighbors can only
//      ever be combined with `partial` alone (it can't join with any other
//      candidate). If v also has no X-neighbor, `partial + {v}` is
//      immediately a maximal clique -- emit it. Either way v is removed
//      from P since it can't participate in any further branching here.
//   2. Lemma 7 (relaxed degree-one in P): a candidate v with exactly one
//      P-neighbor u forms the pair `partial + {v, u}`. This pair is
//      maximal (safe to emit directly) whenever no single x in X is
//      adjacent to *both* v and u -- the code proves this cheaply via the
//      per-call `dynamic_mark_epoch_` marks (see below) instead of
//      recomputing an intersection with X for every such pair.
//   3. Lemma 8, universal candidates: once no more degree-0/1 reductions
//      apply, any candidate adjacent to every *other* remaining candidate
//      must be in every maximal clique of this subproblem, so it's moved
//      from P directly into R (and X is narrowed to its neighborhood, same
//      as a normal BK move-to-R step).
bool Enumerator::apply_rmce_dynamic_reduction(std::vector<int> &partial,
                                              std::vector<int> &candidates,
                                              std::vector<int> &excluded,
                                              int owner_edge) {
  ++result_.counters.dynamic_reduction_calls;
  std::vector<int> removed_candidates;
  if (candidates.empty()) {
    return false;
  }

  // Precompute, for every candidate, whether it has *any* neighbor in X --
  // this is the "marked" bit Lemma 7 needs. Rather than clearing
  // dynamic_mark_epoch_ (sized to the whole graph) on every call, bump a
  // per-call epoch counter and treat "stamped with the current epoch" as
  // the mark; stale marks from earlier calls/epochs just don't match. Only
  // reset the array on the rare uint32_t wraparound.
  //
  // RMCE Algorithm 5 first marks candidates that have any forbidden-set
  // neighbor. Its relaxed degree-one rule fires if at least one endpoint is
  // unmarked, which proves that no x in X extends both endpoints.
  ++dynamic_mark_epoch_value_;
  if (dynamic_mark_epoch_value_ == 0U) {
    std::fill(dynamic_mark_epoch_.begin(), dynamic_mark_epoch_.end(), 0U);
    dynamic_mark_epoch_value_ = 1U;
  }
  for (const int v : candidates) {
    if (intersection_size(excluded, graph_.neighbors(v)) != 0U) {
      dynamic_mark_epoch_[static_cast<std::size_t>(v)] =
          dynamic_mark_epoch_value_;
    }
  }

  // Scan a fixed snapshot of P (scan_order) since the loop body mutates the
  // live `candidates`; the binary_search below skips any v this loop
  // already consumed earlier in the same pass (e.g. as the `u` of an
  // earlier v's degree-one pair).
  const std::vector<int> scan_order = candidates;
  for (const int v : scan_order) {
    if (!std::binary_search(candidates.begin(), candidates.end(), v)) {
      continue;
    }
    const std::vector<int> neighbors_in_p =
        intersect_neighbors(candidates, graph_.neighbors(v));
    if (neighbors_in_p.empty()) {
      // Lemma 5: v has no remaining P-neighbor. Emit partial+{v} only if v
      // also has no X-neighbor (otherwise that X vertex would still extend
      // it, so it wouldn't be maximal); v is removed from P either way.
      if (intersection_size(excluded, graph_.neighbors(v)) == 0U) {
        ++result_.counters.dynamic_degree0_outputs;
        if (options_.collect_cliques || options_.validate_invariants) {
          std::vector<int> clique = partial;
          clique.push_back(v);
          emit_clique(std::move(clique), owner_edge);
        } else {
          record_output_count(partial.size() + 1U);
        }
      }
      erase_sorted(candidates, v);
      removed_candidates.push_back(v);
      ++result_.counters.dynamic_degree0_removed;
      continue;
    }
    if (neighbors_in_p.size() != 1U) {
      continue;
    }

    // Lemma 7: v's sole P-neighbor is u. If both v and u are marked (each
    // has some X-neighbor), we can't yet conclude no single x extends
    // both, so skip -- this pair needs the full BK recursion to resolve.
    // Otherwise, at least one of {v,u} has zero X-neighbors, which means
    // no x in X can be adjacent to both, so partial+{v,u} is maximal.
    const int u = neighbors_in_p[0];
    if (dynamic_mark_epoch_[static_cast<std::size_t>(v)] ==
            dynamic_mark_epoch_value_ &&
        dynamic_mark_epoch_[static_cast<std::size_t>(u)] ==
            dynamic_mark_epoch_value_) {
      continue;
    }
    // Captured before erasing v: u's current P-degree, used just below to
    // decide whether removing v also leaves u with zero P-neighbors (i.e.
    // u's only P-neighbor was v), in which case u is equally spent and can
    // be removed from P too.
    const std::size_t u_degree_before =
        intersection_size(candidates, graph_.neighbors(u));
    ++result_.counters.dynamic_degree1_outputs;
    if (options_.collect_cliques || options_.validate_invariants) {
      std::vector<int> clique = partial;
      clique.push_back(v);
      clique.push_back(u);
      emit_clique(std::move(clique), owner_edge);
    } else {
      record_output_count(partial.size() + 2U);
    }

    if (erase_sorted(candidates, v)) {
      removed_candidates.push_back(v);
      ++result_.counters.dynamic_degree1_removed;
    }
    if (u_degree_before == 1U && erase_sorted(candidates, u)) {
      removed_candidates.push_back(u);
      ++result_.counters.dynamic_degree1_removed;
    }
  }

  // RMCE Lemma 8: a candidate adjacent to every other candidate belongs to
  // every maximal clique of this subproblem and can be moved into R.
  while (!candidates.empty()) {
    int universal = -1;
    for (const int v : candidates) {
      if (intersection_size(candidates, graph_.neighbors(v)) + 1U ==
          candidates.size()) {
        universal = v;
        break;
      }
    }
    if (universal < 0) {
      break;
    }
    partial.push_back(universal);
    erase_sorted(candidates, universal);
    excluded = intersect_neighbors(excluded, graph_.neighbors(universal));
    ++result_.counters.dynamic_universal_moved;
  }

  // A reduced low-degree candidate is deliberately not inserted into X:
  // once a surviving candidate is selected, the low-degree rule proves it
  // cannot extend that branch.  If reduction empties P without selecting a
  // survivor, however, R itself must not be reported when one of the
  // deleted candidates still extends it (e.g., the first K4 edge root).
  for (const int removed : removed_candidates) {
    bool extends_partial = true;
    for (const int v : partial) {
      if (!graph_.adjacent(removed, v)) {
        extends_partial = false;
        break;
      }
    }
    if (extends_partial) {
      return true;
    }
  }
  return false;
}

// RMCE Lemma 9 (forbidden-set containment, Algorithm 6): a vertex x_i in X
// is redundant -- and can be dropped from X without changing the recursion's
// output -- whenever its neighborhood within P, N_P(x_i), is contained in
// some other x_j's N_P(x_j). The reasoning: X's only job is to block
// non-maximal output by proving "some excluded vertex is adjacent to every
// candidate we might select"; if x_i's candidate-neighborhood is a subset
// of x_j's, then any subset of P that x_i would block, x_j also blocks, so
// x_i can never be the sole reason a branch is deemed non-maximal.
//
// When two forbidden vertices have exactly *equal* candidate neighborhoods
// (`equal`), only one is kept -- otherwise both would satisfy `subset` of
// each other and both would be removed, dropping the blocking vertex
// entirely. The `excluded[i] > excluded[j]` tie-break keeps the
// smaller-indexed one deterministically.
void Enumerator::apply_rmce_forbidden_set_reduction(
    const std::vector<int> &candidates, std::vector<int> &excluded) {
  if (excluded.size() < 2U) {
    return;
  }
  // candidate_neighborhoods[i] = N(excluded[i]) ∩ candidates, precomputed
  // once so the O(|X|^2) containment scan below reuses it instead of
  // recomputing an intersection per pair.
  std::vector<std::vector<int>> candidate_neighborhoods;
  candidate_neighborhoods.reserve(excluded.size());
  for (const int x : excluded) {
    candidate_neighborhoods.push_back(
        intersect_neighbors(candidates, graph_.neighbors(x)));
  }

  std::vector<unsigned char> remove(excluded.size(), 0U);
  for (std::size_t i = 0; i < excluded.size(); ++i) {
    for (std::size_t j = 0; j < excluded.size(); ++j) {
      if (i == j) {
        continue;
      }
      const auto &ni = candidate_neighborhoods[i];
      const auto &nj = candidate_neighborhoods[j];
      const bool subset =
          std::includes(nj.begin(), nj.end(), ni.begin(), ni.end());
      const bool equal = ni.size() == nj.size() && subset;
      if (subset && (!equal || excluded[i] > excluded[j])) {
        remove[i] = 1U;
        break;
      }
    }
  }
  std::vector<int> reduced;
  reduced.reserve(excluded.size());
  for (std::size_t i = 0; i < excluded.size(); ++i) {
    if (remove[i] == 0U) {
      reduced.push_back(excluded[i]);
    } else {
      ++result_.counters.forbidden_set_removed;
    }
  }
  excluded = std::move(reduced);
}

// ET only applies at nodes where X is already empty (a prerequisite from
// the paper: with a nonempty X, whether P's cliques are maximal in the
// *original* graph can't be decided from G[P] alone). Tries t=1 up through
// the configured threshold and uses the first (smallest, hence cheapest)
// t for which P is a t-plex; larger t is a strictly weaker condition, so
// this is just picking the tightest applicable case.
bool Enumerator::try_early_termination(const std::vector<int> &partial,
                                       const std::vector<int> &candidates,
                                       const std::vector<int> &excluded,
                                       int owner_edge) {
  if (options_.early_termination_threshold == 0 || !excluded.empty()) {
    return false;
  }
  ++result_.counters.et_checks;

  int accepted_t = 0;
  for (int t = 1; t <= options_.early_termination_threshold; ++t) {
    if (is_t_plex(graph_, candidates, t)) {
      accepted_t = t;
      break;
    }
  }
  if (accepted_t == 0) {
    return false;
  }

  if (accepted_t == 1) {
    ++result_.counters.et1_calls;
  } else if (accepted_t == 2) {
    ++result_.counters.et2_calls;
  } else {
    ++result_.counters.et3_calls;
  }

  // Enumerate every terminal continuation even in count-only mode. The
  // normal output sink still avoids retaining clique identities unless
  // collection or validation was requested, but runtime now includes the
  // per-clique traversal and construction cost instead of an algebraic
  // aggregate count.
  enumerate_t_plex_maximal_cliques(
      graph_, candidates, accepted_t, [&](const std::vector<int> &suffix) {
        std::vector<int> clique = partial;
        clique.insert(clique.end(), suffix.begin(), suffix.end());
        if (accepted_t == 1) {
          ++result_.counters.et1_outputs;
        } else if (accepted_t == 2) {
          ++result_.counters.et2_outputs;
        } else {
          ++result_.counters.et3_outputs;
        }
        emit_clique(std::move(clique), owner_edge);
      });
  return true;
}

// Maximum-candidate-degree pivoting (as named in the README): scans P ∪ X
// and picks whichever vertex has the most neighbors *inside P*, breaking
// ties by smallest vertex id for determinism. Maximizing |N(pivot) ∩ P|
// minimizes |P \ N(pivot)|, i.e. the number of branches vertex_bk's caller
// will have to create, which is the whole point of pivoting.
int Enumerator::choose_pivot(const std::vector<int> &candidates,
                             const std::vector<int> &excluded) const {
  std::vector<int> choices;
  choices.reserve(candidates.size() + excluded.size());
  std::set_union(candidates.begin(), candidates.end(), excluded.begin(),
                 excluded.end(), std::back_inserter(choices));

  int best = -1;
  std::size_t best_neighbors = 0;
  for (const int v : choices) {
    const std::size_t count =
        intersection_size(candidates, graph_.neighbors(v));
    if (best < 0 || count > best_neighbors ||
        (count == best_neighbors && v < best)) {
      best = v;
      best_neighbors = count;
    }
  }
  return best;
}

// True iff (u,v) ranks strictly after `owner_edge` in the truss order --
// the O(1) primitive the entire ownership scheme is built on (see
// enumerate_root_edge and the vertex_bk P/X-split loop). owner_edge==-1
// only occurs for the isolated-vertex case in run(), which has no rank
// context to compare against, so it degrades to plain adjacency.
bool Enumerator::edge_is_after_owner(int u, int v, int owner_edge) const {
  if (owner_edge < 0) {
    return graph_.adjacent(u, v);
  }
  const int edge_id = graph_.edge_id(u, v);
  if (edge_id < 0) {
    return false;
  }
  return result_.truss_order.rank_of_edge[static_cast<std::size_t>(edge_id)] >
         result_.truss_order.rank_of_edge[static_cast<std::size_t>(owner_edge)];
}

// True iff every original-graph edge between two current candidates is
// itself after the owner edge in truss rank -- i.e. G[candidates] (as seen
// by this BK node) is identical to what it would be without any of the
// rank-aware P/X filtering, so RMCE dynamic/forbidden reductions and t-plex
// ET (all designed for an unfiltered induced subgraph) remain sound to
// apply here. This can only be false after at least one vertex_bk
// descent has already demoted an original edge's endpoint from P into X
// (see the "adjacent but earlier-ranked" comment in vertex_bk's branch
// loop); at that point these three optimizations are skipped for the rest
// of the subtree and plain rank-aware BK recursion continues instead (see
// the README's "rank filter" paragraph).
bool Enumerator::candidate_graph_matches_original(
    const std::vector<int> &candidates, int owner_edge) const {
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    for (std::size_t j = i + 1U; j < candidates.size(); ++j) {
      if (graph_.adjacent(candidates[i], candidates[j]) &&
          !edge_is_after_owner(candidates[i], candidates[j], owner_edge)) {
        return false;
      }
    }
  }
  return true;
}

// The count-only fast path: every terminal/ET/dynamic-reduction output goes
// through here (directly, or via emit_clique below) regardless of mode, so
// `maximal_clique_count` is always accurate even when no clique identity is
// ever materialized.
void Enumerator::record_output_count(std::size_t clique_size,
                                     std::uint64_t count) {
  if (clique_size < options_.minimum_output_clique_size) {
    return;
  }
  if (count > std::numeric_limits<std::uint64_t>::max() -
                  result_.maximal_clique_count) {
    throw std::overflow_error("maximal-clique count exceeds uint64_t");
  }
  result_.maximal_clique_count += count;
}

// Identity-materializing output sink. Always tallies the numeric count
// first (record_output_count) regardless of mode; everything below that is
// skipped unless collect_cliques or validate_invariants is set, in which
// case the clique is canonicalized (sorted, deduplicated -- a defensive
// step since ET/dynamic-reduction paths shouldn't produce duplicates, but
// ownership validation below needs a canonical form to look up edges by
// endpoint pairs anyway) and, in validate mode, checked against every
// invariant the README's proof depends on:
//   - is_clique: every pair is actually adjacent in G.
//   - is_maximal: no vertex outside the clique extends it.
//   - clique_key + seen_cliques_: never emitted before (no duplicates
//     across the whole run).
//   - earliest-edge ownership: for size>=2, the clique's own earliest edge
//     (by truss rank) must equal the `owner_edge` this call was made
//     under -- this is the direct empirical check of the ownership
//     invariant the whole algorithm is built on (see enumerator.hpp's
//     class-level comment). Size-1 cliques must have owner_edge==-1 (only
//     the isolated-vertex path in run() is allowed to emit them).
void Enumerator::emit_clique(std::vector<int> clique, int owner_edge) {
  const std::size_t clique_size = clique.size();
  record_output_count(clique_size);
  if (clique_size < options_.minimum_output_clique_size) {
    return;
  }
  if (!options_.collect_cliques && !options_.validate_invariants) {
    return;
  }

  std::sort(clique.begin(), clique.end());
  clique.erase(std::unique(clique.begin(), clique.end()), clique.end());

  if (options_.validate_invariants) {
    if (!is_clique(clique)) {
      ++result_.counters.invalid_clique_outputs;
    }
    if (!is_maximal(clique)) {
      ++result_.counters.nonmaximal_outputs;
    }
    const std::string key = clique_key(clique);
    if (!seen_cliques_.insert(key).second) {
      ++result_.counters.duplicate_outputs;
    }
    if (clique.size() == 1U) {
      if (owner_edge != -1) {
        ++result_.counters.owner_violations;
      }
    } else if (clique.size() >= 2U) {
      // Find the clique's own earliest-ranked internal edge directly from
      // the truss order, independent of how this clique was actually
      // found, and compare it against the owner_edge this call claims.
      int earliest_edge = -1;
      int earliest_rank = std::numeric_limits<int>::max();
      for (std::size_t i = 0; i < clique.size(); ++i) {
        for (std::size_t j = i + 1U; j < clique.size(); ++j) {
          const int edge_id = graph_.edge_id(clique[i], clique[j]);
          if (edge_id < 0) {
            continue;
          }
          const int rank = result_.truss_order
                               .rank_of_edge[static_cast<std::size_t>(edge_id)];
          if (rank < earliest_rank) {
            earliest_rank = rank;
            earliest_edge = edge_id;
          }
        }
      }
      if (earliest_edge != owner_edge) {
        ++result_.counters.owner_violations;
      }
    }
  }

  if (options_.collect_cliques) {
    result_.cliques.push_back(std::move(clique));
  }
}

// Brute-force validation-only checks (only reached under --validate): they
// deliberately don't reuse any state the search built up (P/X, truss rank,
// etc.), so a bug in the fast paths above can't hide a bug from validation.
bool Enumerator::is_clique(const std::vector<int> &clique) const {
  for (std::size_t i = 0; i < clique.size(); ++i) {
    for (std::size_t j = i + 1U; j < clique.size(); ++j) {
      if (!graph_.adjacent(clique[i], clique[j])) {
        return false;
      }
    }
  }
  return true;
}

bool Enumerator::is_maximal(const std::vector<int> &clique) const {
  if (!is_clique(clique)) {
    return false;
  }
  for (int v = 0; v < graph_.vertex_count(); ++v) {
    if (std::binary_search(clique.begin(), clique.end(), v)) {
      continue;
    }
    bool extends = true;
    for (const int u : clique) {
      if (!graph_.adjacent(u, v)) {
        extends = false;
        break;
      }
    }
    if (extends) {
      return false;
    }
  }
  return !clique.empty();
}

std::string Enumerator::clique_key(const std::vector<int> &clique) const {
  std::ostringstream key;
  for (const int v : clique) {
    key << v << ',';
  }
  return key.str();
}

} // namespace hbbmc_faithful
