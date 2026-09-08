#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hbbmc_faithful {

// Undirected edge, always stored/normalized with u < v (see Graph::build).
// `u`/`v` are *dense* vertex ids (0..vertex_count()-1), not original labels.
struct Edge {
    int u = -1;
    int v = -1;
};

// Immutable, simple (no loops/multi-edges) undirected graph.
//
// Every vertex is identified internally by a dense id in [0, vertex_count()),
// which is what the rest of the codebase (truss order, reduction, enumerator)
// operates on. The original, possibly sparse/non-contiguous, input labels
// (e.g. arbitrary integers from a real-world edge-list file) are preserved
// separately in `original_labels_` purely for I/O: they are translated back
// only when printing cliques or looking up input vertices in main.cpp.
//
// Every edge is also assigned a dense edge id in [0, edge_count()), which is
// the unit the truss order and the enumerator's "is this edge after the
// owner edge" checks are built on (see truss_order.hpp, enumerator.cpp).
class Graph {
public:
    Graph() = default;
    // Builds a graph on vertices 0..vertex_count-1 where original labels are
    // just the dense ids themselves (used for synthetic/test graphs).
    Graph(int vertex_count, const std::vector<std::pair<int, int>>& edges);
    // Builds a graph whose dense vertex i corresponds to original_labels[i];
    // `edges` must already be expressed in dense ids. Used when re-labeling
    // a reduced residual graph (see reduction.cpp) so original input labels
    // survive graph reduction.
    Graph(std::vector<std::int64_t> original_labels,
          const std::vector<std::pair<int, int>>& edges);

    // Parses a whitespace-separated edge list ('#'/'%' comment lines and
    // blank lines skipped, self-loops dropped, duplicate edges merged) and
    // assigns dense ids. If declared_vertices >= 0, labels are used directly
    // as dense ids in [0, declared_vertices) so isolated vertices with no
    // incident edge are still retained; otherwise dense ids are assigned only
    // to labels that actually appear in some edge (see read_edge_list impl).
    static Graph read_edge_list(const std::string& path, int declared_vertices = -1);

    int vertex_count() const noexcept { return static_cast<int>(adjacency_.size()); }
    int edge_count() const noexcept { return static_cast<int>(edges_.size()); }

    // Sorted ascending list of dense neighbor ids. Callers throughout the
    // codebase (enumerator, truss order) rely on this being sorted to do
    // set_intersection/set_difference/binary_search on it directly.
    const std::vector<int>& neighbors(int v) const { return adjacency_.at(v); }
    const Edge& edge(int edge_id) const { return edges_.at(edge_id); }
    const std::vector<Edge>& edges() const noexcept { return edges_; }

    bool adjacent(int u, int v) const;
    // Returns the dense edge id of (u,v), or -1 if u==v or the edge is absent.
    int edge_id(int u, int v) const;
    // Sorted ascending vertices adjacent to both u and v (their "support" set
    // used by truss ordering and root-edge candidate/excluded partitioning).
    std::vector<int> common_neighbors(int u, int v) const;

    std::int64_t original_label(int v) const { return original_labels_.at(v); }
    const std::vector<std::int64_t>& original_labels() const noexcept {
        return original_labels_;
    }

private:
    // Order-independent hash key for an unordered pair (u,v): packs the
    // smaller id into the high 32 bits and the larger into the low 32 bits.
    static std::uint64_t edge_key(int u, int v);
    // Normalizes, dedups, and sorts `edges`, assigns dense edge ids in
    // sorted (u,v) order, and builds the adjacency lists. Because edges are
    // processed in sorted order, each adjacency_[x] ends up sorted too (see
    // graph.cpp for the argument), which is what lets adjacent() use
    // binary_search instead of a hash lookup.
    void build(int vertex_count,
               const std::vector<std::pair<int, int>>& edges,
               std::vector<std::int64_t> original_labels);

    std::vector<std::vector<int>> adjacency_;
    std::vector<Edge> edges_;
    std::unordered_map<std::uint64_t, int> edge_ids_;
    std::vector<std::int64_t> original_labels_;
};

}  // namespace hbbmc_faithful
