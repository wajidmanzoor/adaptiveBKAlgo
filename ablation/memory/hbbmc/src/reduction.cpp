#include "../inc/reduction.hpp"

#include <algorithm>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hbbmc_faithful {

namespace {

class NoGraphReduction final : public GraphReductionModule {
public:
  std::string name() const override { return "none"; }
  bool complete_hbbmc_gr() const override { return false; }

  ReductionResult apply(const Graph &graph, ReductionOptions) const override {
    return {graph, 0U, {}, {}};
  }
};

// Implements the RMCE global-reduction pipeline (paper: "Accelerating
// Maximal Clique Enumeration via Graph Reduction", Algorithms 3 and 4):
// repeatedly strip vertices of degree <= 2 and edges that are in no
// triangle, emitting each one's forced maximal clique(s) as they go, until
// no more reductions apply. What's left ("the residual graph") is handed to
// the Enumerator, which only needs to search that smaller graph.
//
// Every rule here only fires when it can *prove* the resulting clique(s)
// are maximal in the *original* graph, which is what makes direct output
// safe to emit without running the BK search at all:
//  - degree-0 vertex present in the original input: an isolated vertex is
//    trivially a maximal 1-clique.
//  - degree-1 vertex u with neighbor v: edge (u,v) can't be extended (u has
//    no other neighbor), so {u,v} is maximal.
//  - degree-2 vertex u with neighbors v,w: if v,w are not adjacent, then
//    {u,v} and {u,w} are both maximal (same argument as degree-1, applied
//    twice); if v,w are adjacent, {u,v,w} is a maximal triangle.
// After emitting, the vertex (and possibly an edge) is deleted, which can
// drop other vertices' degree to <=2 and expose new instances of these
// rules, or make an edge triangle-free -- hence the fixed-point loop below.
class RmceGraphReduction final : public GraphReductionModule {
public:
  std::string name() const override { return "rmce"; }
  bool complete_hbbmc_gr() const override { return true; }

  ReductionResult apply(const Graph &graph,
                        ReductionOptions options) const override {
    const int n = graph.vertex_count();
    // Mutable working copy of the adjacency lists (std::set so individual
    // edge removals during reduction are cheap); `active` distinguishes
    // "vertex fully consumed by a reduction rule" from "vertex still part
    // of the residual graph".
    std::vector<std::set<int>> adjacency(static_cast<std::size_t>(n));
    std::vector<int> original_degree(static_cast<std::size_t>(n), 0);
    std::vector<unsigned char> active(static_cast<std::size_t>(n), 1U);
    for (int v = 0; v < n; ++v) {
      const auto &row = graph.neighbors(v);
      adjacency[static_cast<std::size_t>(v)].insert(row.begin(), row.end());
      original_degree[static_cast<std::size_t>(v)] =
          static_cast<int>(row.size());
    }

    ReductionResult result;
    // Populated only for validation-mode duplicate detection. Every output
    // is retained in directly_emitted_cliques regardless of validation.
    std::set<std::vector<std::int64_t>> emitted;
    auto emit = [&](std::vector<int> dense_clique) {
      if (dense_clique.size() < options.minimum_output_clique_size) {
        return;
      }
      if (result.directly_emitted_count ==
          std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error(
            "direct reduction output count exceeds uint64_t");
      }
      std::vector<std::int64_t> clique;
      clique.reserve(dense_clique.size());
      for (const int v : dense_clique) {
        clique.push_back(graph.original_label(v));
      }
      if (!options.validate_invariants) {
        ++result.directly_emitted_count;
        result.directly_emitted_cliques.push_back(std::move(clique));
        ++result.counters.directly_emitted;
        return;
      }
      std::sort(clique.begin(), clique.end());
      if (emitted.insert(clique).second) {
        ++result.directly_emitted_count;
        result.directly_emitted_cliques.push_back(std::move(clique));
        ++result.counters.directly_emitted;
      } else {
        ++result.counters.duplicate_direct_outputs;
      }
    };

    // Worklist of vertices whose *current* degree is <=2 and therefore due
    // for degree-0/1/2 processing. `queued` prevents a vertex from being
    // pushed twice while it's already waiting in the queue.
    std::queue<int> low_degree;
    std::vector<unsigned char> queued(static_cast<std::size_t>(n), 0U);
    auto enqueue_if_low = [&](int v) {
      if (v >= 0 && active[static_cast<std::size_t>(v)] != 0U &&
          adjacency[static_cast<std::size_t>(v)].size() <= 2U &&
          queued[static_cast<std::size_t>(v)] == 0U) {
        queued[static_cast<std::size_t>(v)] = 1U;
        low_degree.push(v);
      }
    };
    for (int v = 0; v < n; ++v) {
      enqueue_if_low(v);
    }

    // Deletes edge (u,v) from the working adjacency and re-checks both
    // endpoints for the low-degree worklist, since removing an edge can be
    // exactly what drops a neighbor to degree <=2.
    auto remove_edge = [&](int u, int v) {
      adjacency[static_cast<std::size_t>(u)].erase(v);
      adjacency[static_cast<std::size_t>(v)].erase(u);
      enqueue_if_low(u);
      enqueue_if_low(v);
    };

    // Counts common neighbors of u and v among currently-active vertices,
    // via sorted-set merge (std::set iteration is ordered). Used only to
    // decide whether u is v's *sole* remaining common neighbor (see the
    // degree-2 triangle case below).
    auto common_neighbor_count = [&](int u, int v) {
      const auto &a = adjacency[static_cast<std::size_t>(u)];
      const auto &b = adjacency[static_cast<std::size_t>(v)];
      auto i = a.begin();
      auto j = b.begin();
      int count = 0;
      while (i != a.end() && j != b.end()) {
        if (*i < *j) {
          ++i;
        } else if (*j < *i) {
          ++j;
        } else {
          if (active[static_cast<std::size_t>(*i)] != 0U) {
            ++count;
          }
          ++i;
          ++j;
        }
      }
      return count;
    };

    // Drains the low-degree worklist, applying the degree-0/1/2 vertex
    // rules (see the class-level comment above) until every currently
    // low-degree vertex has been consumed. Returns whether anything changed
    // this call, which the outer fixed-point loop uses to decide whether
    // another pass over triangle-free edges is worth trying.
    auto process_low_degree = [&]() {
      bool changed = false;
      while (!low_degree.empty()) {
        const int u = low_degree.front();
        low_degree.pop();
        queued[static_cast<std::size_t>(u)] = 0U;
        if (active[static_cast<std::size_t>(u)] == 0U ||
            adjacency[static_cast<std::size_t>(u)].size() > 2U) {
          continue;
        }
        changed = true;
        const std::vector<int> neighbors(
            adjacency[static_cast<std::size_t>(u)].begin(),
            adjacency[static_cast<std::size_t>(u)].end());

        if (neighbors.empty()) {
          ++result.counters.degree0_vertices;
          // RMCE's degree-zero lemma assumes cliques have size >=2.
          // HBBMC explicitly includes maximal 1-cliques, so an
          // isolate present in the original graph is emitted here.
          // Vertices made isolated by earlier reductions are silent:
          // their incident maximal cliques were already reported.
          if (original_degree[static_cast<std::size_t>(u)] == 0) {
            emit({u});
          }
          active[static_cast<std::size_t>(u)] = 0U;
          continue;
        }

        if (neighbors.size() == 1U) {
          ++result.counters.degree1_vertices;
          const int v = neighbors[0];
          emit({u, v});
          remove_edge(u, v);
          active[static_cast<std::size_t>(u)] = 0U;
          adjacency[static_cast<std::size_t>(u)].clear();
          enqueue_if_low(v);
          continue;
        }

        ++result.counters.degree2_vertices;
        const int v = neighbors[0];
        const int w = neighbors[1];
        const bool vw_edge = adjacency[static_cast<std::size_t>(v)].find(w) !=
                             adjacency[static_cast<std::size_t>(v)].end();
        if (!vw_edge) {
          emit({u, v});
          emit({u, w});
        } else {
          emit({u, v, w});
          // RMCE Lemma 3, case 2: when u is the only common
          // neighbor, deleting (v,w) prevents a duplicate 2-clique.
          if (common_neighbor_count(v, w) == 1) {
            remove_edge(v, w);
          }
        }
        remove_edge(u, v);
        remove_edge(u, w);
        active[static_cast<std::size_t>(u)] = 0U;
        adjacency[static_cast<std::size_t>(u)].clear();
        enqueue_if_low(v);
        enqueue_if_low(w);
      }
      return changed;
    };

    // RMCE Algorithms 3 and 4. We repeat the two phases because deleting
    // non-triangle edges can expose new degree-at-most-two vertices, as in
    // the paper's Figure 3 example.
    while (true) {
      // Phase 1 (Algorithm 3): drain every currently low-degree vertex.
      bool changed = process_low_degree();

      // Phase 2 (Algorithm 4): find every remaining edge with zero common
      // neighbors (i.e. in no triangle) -- such an edge cannot be extended
      // by any vertex, so it is itself a maximal 2-clique. Collect the
      // candidate edge list first (snapshotting u<v pairs) since remove_edge
      // mutates `adjacency` and would otherwise invalidate iteration.
      std::vector<std::pair<int, int>> edges;
      for (int u = 0; u < n; ++u) {
        if (active[static_cast<std::size_t>(u)] == 0U) {
          continue;
        }
        for (const int v : adjacency[static_cast<std::size_t>(u)]) {
          if (u < v && active[static_cast<std::size_t>(v)] != 0U) {
            edges.emplace_back(u, v);
          }
        }
      }
      for (const auto &[u, v] : edges) {
        // Re-check liveness: an earlier iteration of this same loop may
        // already have deleted (u,v) or one of its endpoints.
        if (active[static_cast<std::size_t>(u)] == 0U ||
            active[static_cast<std::size_t>(v)] == 0U ||
            adjacency[static_cast<std::size_t>(u)].find(v) ==
                adjacency[static_cast<std::size_t>(u)].end()) {
          continue;
        }
        if (common_neighbor_count(u, v) == 0) {
          emit({u, v});
          ++result.counters.nontriangle_edges;
          remove_edge(u, v);
          changed = true;
        }
      }
      // Fixed point: stop only once a full pass produced no new low-degree
      // work and no triangle-free edge was found.
      if (!changed && low_degree.empty()) {
        break;
      }
    }

    // Build the residual graph from vertices that are still active. The
    // low-degree worklist is fully drained on every process_low_degree()
    // call (including entries it enqueues while running), so by the time
    // the fixed-point loop above exits, every active vertex should already
    // have degree > 2; the `!adjacency.empty()` guard is a defensive
    // belt-and-suspenders check against emitting a stray isolated vertex
    // into the residual graph rather than a reachable steady-state case.
    std::vector<int> old_to_new(static_cast<std::size_t>(n), -1);
    std::vector<std::int64_t> labels;
    for (int old = 0; old < n; ++old) {
      if (active[static_cast<std::size_t>(old)] != 0U &&
          !adjacency[static_cast<std::size_t>(old)].empty()) {
        old_to_new[static_cast<std::size_t>(old)] =
            static_cast<int>(labels.size());
        labels.push_back(graph.original_label(old));
      }
    }
    std::vector<std::pair<int, int>> residual_edges;
    for (int old_u = 0; old_u < n; ++old_u) {
      const int u = old_to_new[static_cast<std::size_t>(old_u)];
      if (u < 0) {
        continue;
      }
      for (const int old_v : adjacency[static_cast<std::size_t>(old_u)]) {
        const int v = old_to_new[static_cast<std::size_t>(old_v)];
        if (v >= 0 && u < v) {
          residual_edges.emplace_back(u, v);
        }
      }
    }
    result.graph = Graph(std::move(labels), residual_edges);
    return result;
  }
};

} // namespace

std::unique_ptr<GraphReductionModule>
make_graph_reduction_module(const std::string &name) {
  if (name == "none") {
    return std::make_unique<NoGraphReduction>();
  }
  if (name == "rmce" || name == "gr") {
    return std::make_unique<RmceGraphReduction>();
  }
  throw std::invalid_argument("unknown graph-reduction module: " + name);
}

} // namespace hbbmc_faithful
