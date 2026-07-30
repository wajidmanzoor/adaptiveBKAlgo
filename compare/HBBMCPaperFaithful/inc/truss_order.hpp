#pragma once

#include "graph.hpp"

#include <vector>

namespace hbbmc_faithful {

struct TrussOrder {
  // edge_at_rank[r] is the edge removed at rank r (zero based).
  std::vector<int> edge_at_rank;
  std::vector<int> rank_of_edge;
  std::vector<int> support_at_removal;
  int tau = 0;
};

// Implements the paper's greedy truss-based ordering: repeatedly remove a
// remaining edge with minimum current support (number of common neighbors).
// Ties are resolved by the stable edge id for deterministic experiments.
TrussOrder compute_truss_order(const Graph &graph);

} // namespace hbbmc_faithful
