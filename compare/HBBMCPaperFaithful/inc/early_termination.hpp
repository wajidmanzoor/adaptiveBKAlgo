#pragma once

#include "graph.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace hbbmc_faithful {

using CandidateCliqueCallback = std::function<void(const std::vector<int> &)>;

// A candidate graph is a t-plex under the convention used in the HBBMC
// paper when every candidate has at most t-1 non-neighbors in candidates.
bool is_t_plex(const Graph &graph, const std::vector<int> &candidates, int t);

// Algebraic cross-check helper. The production Enumerator intentionally uses
// enumerate_t_plex_maximal_cliques even in count-only mode so its runtime
// includes one visit and construction per terminal clique. Throws on uint64_t
// overflow.
std::uint64_t count_t_plex_maximal_cliques(const Graph &graph,
                                           const std::vector<int> &candidates,
                                           int t);

// Enumerate every maximal clique of G[candidates] without BK branching.
// Preconditions: t is 1, 2, or 3 and candidates induce a t-plex.
void enumerate_t_plex_maximal_cliques(const Graph &graph,
                                      const std::vector<int> &candidates, int t,
                                      const CandidateCliqueCallback &callback);

} // namespace hbbmc_faithful
