#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hbbmc_faithful {

struct Edge {
    int u = -1;
    int v = -1;
};

class Graph {
public:
    Graph() = default;
    Graph(int vertex_count, const std::vector<std::pair<int, int>>& edges);
    Graph(std::vector<std::int64_t> original_labels,
          const std::vector<std::pair<int, int>>& edges);

    static Graph read_edge_list(const std::string& path, int declared_vertices = -1);

    int vertex_count() const noexcept { return static_cast<int>(adjacency_.size()); }
    int edge_count() const noexcept { return static_cast<int>(edges_.size()); }

    const std::vector<int>& neighbors(int v) const { return adjacency_.at(v); }
    const Edge& edge(int edge_id) const { return edges_.at(edge_id); }
    const std::vector<Edge>& edges() const noexcept { return edges_; }

    bool adjacent(int u, int v) const;
    int edge_id(int u, int v) const;
    std::vector<int> common_neighbors(int u, int v) const;

    std::int64_t original_label(int v) const { return original_labels_.at(v); }
    const std::vector<std::int64_t>& original_labels() const noexcept {
        return original_labels_;
    }

private:
    static std::uint64_t edge_key(int u, int v);
    void build(int vertex_count,
               const std::vector<std::pair<int, int>>& edges,
               std::vector<std::int64_t> original_labels);

    std::vector<std::vector<int>> adjacency_;
    std::vector<Edge> edges_;
    std::unordered_map<std::uint64_t, int> edge_ids_;
    std::vector<std::int64_t> original_labels_;
};

}  // namespace hbbmc_faithful
