#pragma once

#include "fast_adj_hash.h"
#include "fast_clique_sink.h"

struct FastLocalBitsetResult {
  bool handled = false;
  bool found = false;
  ull cliqueCount = 0;
  ui maxCliqueSize = 0;
  ull checksCount = 0;
  std::vector<ui> witness;
};

// Solve one Bron--Kerbosch subtree in a single machine word. The input state
// is read-only: P and X keep their full BK meaning, cliqueSize is |R|, and
// cliquePrefix is the materialized R used only when a sibling witness is
// requested. When cliqueSink is supplied, cliquePrefix must contain R and
// every output is sent to the sink in canonical sorted order. handled=false
// asks the caller to retain its list recursion.
FastLocalBitsetResult solveFastLocalBitsetSubtree(
    const FastAdjacencyHash &adjacency, const std::vector<ui> &p,
    const std::vector<ui> &x, ui cliqueSize,
    const std::vector<ui> *cliquePrefix, ui minCliqueSize = 3,
    const FastCliqueSink *cliqueSink = nullptr);
