#pragma once

#include "fast_clique_sink.h"

#include <unordered_set>

struct FastPlex3Result {
  bool handled = false;
  bool found = false;
  ull cliqueCount = 0;
  ui maxCliqueSize = 0;
  ull checksCount = 0;
  std::vector<ui> witness;
};

// Solve an X-empty Bron--Kerbosch subtree when the complement of P has
// maximum degree at most two (equivalently, G[P] is a 3-plex).  Maximal
// clique continuations are maximal independent sets in that complement,
// whose connected components are isolated vertices, paths, or cycles.
//
// cliqueSize is |R|. When cliquePrefix is supplied, it is prepended to each
// materialized clique. handled=false asks the caller to retain ordinary
// recursion. It is returned for non-3-plex states and whenever a result would
// overflow the public counters.
FastPlex3Result solveFastPlex3Subtree(
    const std::vector<std::unordered_set<ui>> &adjacency,
    const std::vector<ui> &p, ui cliqueSize,
    const std::vector<ui> *cliquePrefix = nullptr, ui minCliqueSize = 3,
    const FastCliqueSink *cliqueSink = nullptr);
