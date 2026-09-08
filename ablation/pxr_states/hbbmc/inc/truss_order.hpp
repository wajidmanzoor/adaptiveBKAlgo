#pragma once

#include "graph.hpp"

#include <vector>

namespace hbbmc_faithful {

// The output of truss peeling: a total order over all edges, produced by
// repeatedly deleting the currently-lowest-support edge (support = number of
// triangles the edge participates in, i.e. |common_neighbors(u,v)|).
//
// This order is the backbone of the whole algorithm: every maximal clique of
// size >= 2 is later assigned to its *earliest* edge in this order (its
// "owner"), and the enumerator only ever needs to search forward from an
// edge's rank (see Enumerator::enumerate_root_edge / edge_is_after_owner).
// `tau` (the max support seen at any removal) is the graph's truss number
// and upper-bounds how large a root's candidate set can be.
struct TrussOrder {
  // edge_at_rank[r] is the dense edge id removed at rank r (zero based); the
  // enumerator iterates roots in this order (rank 0 first).
  std::vector<int> edge_at_rank;
  // Inverse of edge_at_rank: rank_of_edge[edge_id] is that edge's removal
  // rank. Used everywhere to compare "is edge A after edge B" in O(1).
  std::vector<int> rank_of_edge;
  // support_at_removal[edge_id] is the number of common neighbors the edge
  // had at the moment it was removed (i.e. its k-truss support), used only
  // for a diagnostic cross-check against the root candidate-set size.
  std::vector<int> support_at_removal;
  // Maximum support_at_removal across all edges; the graph's truss number.
  int tau = 0;
};

// Implements the paper's greedy truss-based ordering: repeatedly remove a
// remaining edge with minimum current support (number of common neighbors).
// Ties are resolved by the stable edge id for deterministic experiments.
TrussOrder compute_truss_order(const Graph &graph);

} // namespace hbbmc_faithful
