#pragma once

#include "graph.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace hbbmc_faithful {

// Why "early termination" exists at all: once BK recursion reaches a node
// where X (excluded) is empty and P (candidates) induces a t-plex (defined
// below) for some t<=3, the paper shows the recursion can stop branching
// vertex-by-vertex and instead directly enumerate every maximal clique of
// G[P] using closed-form recurrences on the *complement* of G[P]. This is
// cheaper than continuing pivoted BK because a t-plex's complement graph has
// very restricted structure (see enumerate_t_plex_maximal_cliques below),
// which collapses what would be branching search into direct construction.
using CandidateCliqueCallback = std::function<void(const std::vector<int> &)>;

// A candidate graph is a t-plex under the convention used in the HBBMC
// paper when every candidate has at most t-1 non-neighbors in candidates.
// Equivalently: every vertex in the complement of G[candidates] has degree
// <= t-1. t=1 means G[candidates] is itself a complete graph (a clique).
bool is_t_plex(const Graph &graph, const std::vector<int> &candidates, int t);

// Algebraic cross-check helper: computes the same count that
// enumerate_t_plex_maximal_cliques would produce, but via closed-form
// component-size formulas (see path_option_count/cycle_option_count in the
// .cpp) instead of actually constructing each clique. This exists only to
// validate the enumeration path in tests; the production Enumerator always
// uses enumerate_t_plex_maximal_cliques so it visits, constructs, and retains
// every terminal clique. Throws on uint64_t overflow.
std::uint64_t count_t_plex_maximal_cliques(const Graph &graph,
                                           const std::vector<int> &candidates,
                                           int t);

// Enumerate every maximal clique of G[candidates] without BK branching.
// Preconditions: t is 1, 2, or 3 and candidates induce a t-plex.
//
// Core idea: a maximal clique of G[candidates] is exactly the complement of
// a maximal *independent set* of the complement graph \overline{G[candidates]}.
// Because that complement graph has max degree t-1 <= 2, it decomposes into
// disjoint connected components that are each a single vertex, a path, or
// (only possible at t=3) a cycle -- see build_complement/is_cycle in the
// .cpp. Each component's maximal independent sets can be enumerated (or
// counted) independently via small direct-recursion/closed-form rules, and
// the results are combined across components by a cartesian product (see
// enumerate_product), since choices in one component don't affect another.
void enumerate_t_plex_maximal_cliques(const Graph &graph,
                                      const std::vector<int> &candidates, int t,
                                      const CandidateCliqueCallback &callback);

} // namespace hbbmc_faithful
