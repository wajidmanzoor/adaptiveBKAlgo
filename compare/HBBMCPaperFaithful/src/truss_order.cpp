#include "inc/truss_order.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace hbbmc_faithful {

TrussOrder compute_truss_order(const Graph &graph) {
  const int n = graph.vertex_count();
  const int m = graph.edge_count();

  TrussOrder result;
  result.edge_at_rank.reserve(static_cast<std::size_t>(m));
  result.rank_of_edge.assign(static_cast<std::size_t>(m), -1);
  result.support_at_removal.assign(static_cast<std::size_t>(m), 0);
  if (m == 0) {
    return result;
  }

  std::vector<std::unordered_set<int>> alive_neighbors(
      static_cast<std::size_t>(n));
  for (int v = 0; v < n; ++v) {
    const auto &row = graph.neighbors(v);
    alive_neighbors[static_cast<std::size_t>(v)].insert(row.begin(), row.end());
  }

  std::vector<int> support(static_cast<std::size_t>(m), 0);
  std::vector<unsigned char> alive(static_cast<std::size_t>(m), 1U);
  using HeapItem = std::pair<int, int>;
  std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>>
      heap;
  for (int edge_id = 0; edge_id < m; ++edge_id) {
    const auto &edge = graph.edge(edge_id);
    support[static_cast<std::size_t>(edge_id)] =
        static_cast<int>(graph.common_neighbors(edge.u, edge.v).size());
    heap.emplace(support[static_cast<std::size_t>(edge_id)], edge_id);
  }

  while (static_cast<int>(result.edge_at_rank.size()) < m) {
    while (!heap.empty()) {
      const auto [value, edge_id] = heap.top();
      if (alive[static_cast<std::size_t>(edge_id)] != 0U &&
          support[static_cast<std::size_t>(edge_id)] == value) {
        break;
      }
      heap.pop();
    }
    if (heap.empty()) {
      throw std::logic_error("truss-order heap exhausted before all edges");
    }

    const auto [value, edge_id] = heap.top();
    heap.pop();
    const Edge edge = graph.edge(edge_id);
    const int rank = static_cast<int>(result.edge_at_rank.size());
    result.edge_at_rank.push_back(edge_id);
    result.rank_of_edge[static_cast<std::size_t>(edge_id)] = rank;
    result.support_at_removal[static_cast<std::size_t>(edge_id)] = value;
    result.tau = std::max(result.tau, value);

    auto &nu = alive_neighbors[static_cast<std::size_t>(edge.u)];
    auto &nv = alive_neighbors[static_cast<std::size_t>(edge.v)];
    const auto *smaller = &nu;
    const auto *larger = &nv;
    if (smaller->size() > larger->size()) {
      std::swap(smaller, larger);
    }
    for (const int w : *smaller) {
      if (w == edge.u || w == edge.v || larger->find(w) == larger->end()) {
        continue;
      }
      const int uw = graph.edge_id(edge.u, w);
      const int vw = graph.edge_id(edge.v, w);
      if (uw < 0 || vw < 0 || alive[static_cast<std::size_t>(uw)] == 0U ||
          alive[static_cast<std::size_t>(vw)] == 0U) {
        continue;
      }
      --support[static_cast<std::size_t>(uw)];
      --support[static_cast<std::size_t>(vw)];
      if (support[static_cast<std::size_t>(uw)] < 0 ||
          support[static_cast<std::size_t>(vw)] < 0) {
        throw std::logic_error("negative support during truss peeling");
      }
      heap.emplace(support[static_cast<std::size_t>(uw)], uw);
      heap.emplace(support[static_cast<std::size_t>(vw)], vw);
    }

    alive[static_cast<std::size_t>(edge_id)] = 0U;
    nu.erase(edge.v);
    nv.erase(edge.u);
  }

  return result;
}

} // namespace hbbmc_faithful
