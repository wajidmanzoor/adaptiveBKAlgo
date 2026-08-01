#include "../inc/enumerator.hpp"

#include "../inc/early_termination.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace hbbmc_faithful {

namespace {

std::vector<int> intersect_neighbors(const std::vector<int> &vertices,
                                     const std::vector<int> &neighbors) {
  std::vector<int> output;
  output.reserve(std::min(vertices.size(), neighbors.size()));
  std::set_intersection(vertices.begin(), vertices.end(), neighbors.begin(),
                        neighbors.end(), std::back_inserter(output));
  return output;
}

std::vector<int> subtract_neighbors(const std::vector<int> &vertices,
                                    const std::vector<int> &neighbors) {
  std::vector<int> output;
  output.reserve(vertices.size());
  std::set_difference(vertices.begin(), vertices.end(), neighbors.begin(),
                      neighbors.end(), std::back_inserter(output));
  return output;
}

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

void insert_sorted_unique(std::vector<int> &values, int value) {
  const auto position = std::lower_bound(values.begin(), values.end(), value);
  if (position == values.end() || *position != value) {
    values.insert(position, value);
  }
}

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
  result_.truss_order = compute_truss_order(graph_);

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

  for (std::size_t rank = 0; rank < result_.truss_order.edge_at_rank.size();
       ++rank) {
    enumerate_root_edge(result_.truss_order.edge_at_rank[rank],
                        static_cast<int>(rank));
  }
  return std::move(result_);
}

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

void Enumerator::vertex_bk(std::vector<int> &partial,
                           std::vector<int> candidates,
                           std::vector<int> excluded, int owner_edge,
                           std::uint64_t depth) {
  const std::size_t partial_size_on_entry = partial.size();
  struct PartialRestore {
    std::vector<int> &partial;
    std::size_t size;
    ~PartialRestore() { partial.resize(size); }
  } partial_restore{partial, partial_size_on_entry};

  ++result_.counters.vertex_recursive_calls;
  result_.counters.max_vertex_depth =
      std::max(result_.counters.max_vertex_depth, depth);

  // RMCE applies maximality-check reduction before its recursive kernel and
  // dynamic reduction at kernel entry. In HBBMC, each truss-owned edge child
  // is the analogous top-level subproblem.
  const bool owner_graph_unfiltered =
      candidate_graph_matches_original(candidates, owner_edge);
  if (options_.rmce_forbidden_set_reduction && depth == 0U &&
      owner_graph_unfiltered) {
    apply_rmce_forbidden_set_reduction(candidates, excluded);
  }
  bool partial_blocked_by_reduced_candidate = false;
  if (options_.rmce_dynamic_reduction && owner_graph_unfiltered) {
    partial_blocked_by_reduced_candidate =
        apply_rmce_dynamic_reduction(partial, candidates, excluded, owner_edge);
  }

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

  if (owner_graph_unfiltered &&
      try_early_termination(partial, candidates, excluded, owner_edge)) {
    return;
  }

  const int pivot = choose_pivot(candidates, excluded);
  if (pivot < 0) {
    throw std::logic_error("nonempty candidate set has no pivot");
  }
  ++result_.counters.pivot_calls;
  const std::vector<int> branch_vertices =
      subtract_neighbors(candidates, graph_.neighbors(pivot));

  for (const int v : branch_vertices) {
    const auto position =
        std::lower_bound(candidates.begin(), candidates.end(), v);
    if (position == candidates.end() || *position != v) {
      continue;
    }
    ++result_.counters.pivot_branches;
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

    candidates.erase(std::lower_bound(candidates.begin(), candidates.end(), v));
    insert_sorted_unique(excluded, v);
  }
}

bool Enumerator::apply_rmce_dynamic_reduction(std::vector<int> &partial,
                                              std::vector<int> &candidates,
                                              std::vector<int> &excluded,
                                              int owner_edge) {
  ++result_.counters.dynamic_reduction_calls;
  std::vector<int> removed_candidates;
  if (candidates.empty()) {
    return false;
  }

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

  const std::vector<int> scan_order = candidates;
  for (const int v : scan_order) {
    if (!std::binary_search(candidates.begin(), candidates.end(), v)) {
      continue;
    }
    const std::vector<int> neighbors_in_p =
        intersect_neighbors(candidates, graph_.neighbors(v));
    if (neighbors_in_p.empty()) {
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

    const int u = neighbors_in_p[0];
    if (dynamic_mark_epoch_[static_cast<std::size_t>(v)] ==
            dynamic_mark_epoch_value_ &&
        dynamic_mark_epoch_[static_cast<std::size_t>(u)] ==
            dynamic_mark_epoch_value_) {
      continue;
    }
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

void Enumerator::apply_rmce_forbidden_set_reduction(
    const std::vector<int> &candidates, std::vector<int> &excluded) {
  if (excluded.size() < 2U) {
    return;
  }
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
