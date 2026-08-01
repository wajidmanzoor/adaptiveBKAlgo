#include "../inc/fast_list_bk.h"
#include "../inc/fast_local_bitset.h"
#include "../inc/fast_plex3.h"

#include <chrono>
#include <iomanip>

#ifdef FASTLIST_OPPORTUNITY_PROFILE
namespace {
ull profileMix(ull value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}
} // namespace

size_t FastListBK::profileBucket(ull value) {
  if (value <= 4)
    return static_cast<size_t>(value);
  if (value <= 8)
    return 5;
  if (value <= 16)
    return 6;
  if (value <= 32)
    return 7;
  if (value <= 64)
    return 8;
  if (value <= 128)
    return 9;
  if (value <= 256)
    return 10;
  if (value <= 512)
    return 11;
  return 12;
}

ull FastListBK::profilePivotItemsTotal() const {
  return profile.pivotHashItemsX + profile.pivotHashItemsP +
         profile.pivotCsrItemsX + profile.pivotCsrItemsP;
}

ull FastListBK::profileIntersectItemsTotal() const {
  return profile.intersectHashItems + profile.intersectCsrItems;
}

void FastListBK::profileGraphRules() {
  for (ui degree : graph.degree) {
    profile.graphDegree0 += degree == 0;
    profile.graphDegree1 += degree == 1;
    profile.graphDegree2 += degree == 2;
  }

  for (ui u = 0; u < graph.n; ++u) {
    for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at) {
      const ui v = graph.neighbors[at];
      if (u >= v)
        continue;
      ++profile.graphUndirectedEdges;

      const ui scan = graph.degree[u] <= graph.degree[v] ? u : v;
      const ui probe = scan == u ? v : u;
      bool inTriangle = false;
      for (ui pos = graph.offset[scan]; pos < graph.offset[scan + 1]; ++pos) {
        const ui w = graph.neighbors[pos];
        if (w == probe)
          continue;
        ++profile.graphNonTriangleProbes;
        if (adjacency.contains(probe, w)) {
          inTriangle = true;
          break;
        }
      }
      profile.graphNonTriangleEdges += !inTriangle;
    }
  }
}

void FastListBK::profileXDominance(ui depth, const Level &level) {
  const size_t pSize = level.p.size();
  const size_t xSize = level.x.size();
  if (pSize == 0 || xSize < 2 || pSize + xSize > 64)
    return;

  ++profile.xDomEligibleSmallNodes;
  const bool warmup = profile.xDomEligibleSmallNodes <= 64;
  const ull sampleKey =
      profileMix(checksCount ^ (static_cast<ull>(depth) << 48));
  if (!warmup && (sampleKey & 8191ULL) != 0)
    return;

  ++profile.xDomSampledNodes;
  profile.xDomWarmupNodes += warmup;
  profile.xDomBefore += xSize;

  std::vector<ull> masks(xSize, 0);
  for (size_t xi = 0; xi < xSize; ++xi) {
    ull mask = 0;
    for (size_t pi = 0; pi < pSize; ++pi) {
      if (adjacency.contains(level.x[xi], level.p[pi]))
        mask |= 1ULL << pi;
    }
    masks[xi] = mask;
  }

  ull removed = 0;
  for (size_t i = 0; i < xSize; ++i) {
    bool dominated = false;
    const int popI = __builtin_popcountll(masks[i]);
    for (size_t j = 0; j < xSize; ++j) {
      if (i == j || __builtin_popcountll(masks[j]) < popI)
        continue;
      ++profile.xDomSubsetTests;
      if ((masks[i] & ~masks[j]) == 0 &&
          (masks[i] != masks[j] || j < i)) {
        dominated = true;
        break;
      }
    }
    removed += dominated;
  }

  profile.xDomRemoved += removed;
  profile.xDomNodesWithRemoval += removed != 0;
}

void FastListBK::printOpportunityProfile() const {
  std::cout << "OPROF_META"
            << " nodes=" << checksCount
            << " degeneracy=" << degeneracy
            << " sibling=" << hybridReorderSibling << std::endl;
  std::cout << "OPROF_GLOBAL"
            << " degree0=" << profile.graphDegree0
            << " degree1=" << profile.graphDegree1
            << " degree2=" << profile.graphDegree2
            << " edges=" << profile.graphUndirectedEdges
            << " nontriangle_edges=" << profile.graphNonTriangleEdges
            << " nontriangle_probes=" << profile.graphNonTriangleProbes
            << std::endl;
  std::cout << "OPROF_TERMINALS"
            << " leaf_maximal=" << profile.leafMaximal
            << " leaf_blocked_x=" << profile.leafBlockedX
            << " x_universal=" << profile.xUniversalPrune
            << " p_clique=" << profile.pCliqueSolved
            << " matching=" << profile.matchingSolved
            << " matching_overflow=" << profile.matchingOverflowFallback
            << " tiny_nodes=" << profile.tinyKernelNodes
            << " tiny_cliques=" << profile.tinyKernelCliques
            << " ordinary_nodes=" << profile.ordinaryBranchNodes
            << " ordinary_branches=" << profile.ordinaryBranchVertices
            << " ordinary_tiny3=" << profile.ordinaryTiny3Nodes << std::endl;
  std::cout << "OPROF_RULES"
            << " d0_nodes=" << profile.pDeg0Nodes
            << " d0_vertices=" << profile.pDeg0Vertices
            << " d0_outside_tiny3=" << profile.pDeg0OutsideTiny3
            << " d1_nodes=" << profile.pDeg1Nodes
            << " d1_vertices=" << profile.pDeg1Vertices
            << " d1_outside_tiny3=" << profile.pDeg1OutsideTiny3
            << " universal_nodes=" << profile.universalPNodes
            << " universal_vertices=" << profile.universalPVertices
            << " universal_outside_tiny3="
            << profile.universalPOutsideTiny3 << std::endl;
  std::cout << "OPROF_PLEX3"
            << " nodes=" << profile.plex3Nodes
            << " p_vertices=" << profile.plex3PVertices
            << " complement_edges=" << profile.plex3ComplementEdges
            << " topmost_nodes=" << profile.plex3TopmostNodes
            << " descendant_checks=" << profile.plex3DescendantChecks
            << " descendant_cliques=" << profile.plex3DescendantCliques
            << " descendant_pivot_items="
            << profile.plex3DescendantPivotItems
            << " subtree_intersect_items="
            << profile.plex3SubtreeIntersectItems << std::endl;
  std::cout << "OPROF_PIVOT"
            << " candidates_x=" << profile.pivotCandidatesX
            << " candidates_p=" << profile.pivotCandidatesP
            << " hash_items_x=" << profile.pivotHashItemsX
            << " hash_items_p=" << profile.pivotHashItemsP
            << " csr_items_x=" << profile.pivotCsrItemsX
            << " csr_items_p=" << profile.pivotCsrItemsP
            << " abortable_x=" << profile.pivotAbortableX
            << " abortable_p=" << profile.pivotAbortableP
            << " abort_saved_x=" << profile.pivotAbortSavedX
            << " abort_saved_p=" << profile.pivotAbortSavedP
            << " naude_width2_nodes=" << profile.naudeWidth2Nodes
            << " naude_tail_candidates=" << profile.naudeTailCandidates
            << " naude_tail_items=" << profile.naudeTailItems
            << " p_improved_nodes=" << profile.pImprovedPivotNodes
            << " winner_x=" << profile.pivotWinnerX
            << " winner_p=" << profile.pivotWinnerP << std::endl;
  std::cout << "OPROF_INTERSECT"
            << " hash_calls=" << profile.intersectHashCalls
            << " hash_items=" << profile.intersectHashItems
            << " csr_calls=" << profile.intersectCsrCalls
            << " csr_items=" << profile.intersectCsrItems << std::endl;
  std::cout << "OPROF_XDOM"
            << " eligible_small_nodes=" << profile.xDomEligibleSmallNodes
            << " sampled_nodes=" << profile.xDomSampledNodes
            << " warmup_nodes=" << profile.xDomWarmupNodes
            << " nodes_with_removal=" << profile.xDomNodesWithRemoval
            << " x_before=" << profile.xDomBefore
            << " x_removed=" << profile.xDomRemoved
            << " subset_tests=" << profile.xDomSubsetTests
            << " random_sample_rate=8192"
            << " state_limit=64" << std::endl;

  for (size_t i = 0; i < PROFILE_BUCKET_COUNT; ++i) {
    std::cout << "OPROF_BUCKET"
              << " bucket=" << i << " nodes=" << profile.stateNodes[i]
              << " pivot_items=" << profile.pivotItems[i]
              << " intersect_items=" << profile.intersectItems[i]
              << std::endl;
  }
  for (size_t i = 0; i < PROFILE_THRESHOLD_COUNT; ++i) {
    const ull fullCovered =
        profile.fullEligible[i] - profile.fullEmptyCrossings[i];
    const ull partialCovered =
        profile.partialEligible[i] - profile.partialEmptyCrossings[i];
    std::cout << "OPROF_CUTOVER"
              << " kind=full threshold=" << PROFILE_THRESHOLDS[i]
              << " eligible=" << profile.fullEligible[i]
              << " empty_crossings=" << profile.fullEmptyCrossings[i]
              << " frontiers=" << profile.fullFrontiers[i]
              << " descendant_calls="
              << (fullCovered - profile.fullFrontiers[i])
              << " frontier_p=" << profile.fullFrontierP[i]
              << " frontier_x=" << profile.fullFrontierX[i]
              << " cells=" << profile.fullCells[i]
              << " words=" << profile.fullWords[i] << std::endl;
    std::cout << "OPROF_CUTOVER"
              << " kind=partial threshold=" << PROFILE_THRESHOLDS[i]
              << " eligible=" << profile.partialEligible[i]
              << " empty_crossings=" << profile.partialEmptyCrossings[i]
              << " frontiers=" << profile.partialFrontiers[i]
              << " descendant_calls="
              << (partialCovered - profile.partialFrontiers[i])
              << " frontier_p=" << profile.partialFrontierP[i]
              << " frontier_x=" << profile.partialFrontierX[i]
              << " cells=" << profile.partialCells[i]
              << " words=" << profile.partialWords[i] << std::endl;
  }
  for (size_t p = 0; p <= 16; ++p) {
    for (size_t x = 0; x <= 16; ++x) {
      const ull nodes = profile.tinyPX[p * 17 + x];
      if (nodes != 0)
        std::cout << "OPROF_TINYPX"
                  << " p=" << p << " x=" << x << " nodes=" << nodes
                  << std::endl;
    }
  }
}
#endif

FastListBK::FastListBK(const Graph &g, bool useHybridReorderSibling,
                       ui outputThreshold)
    : graph(g), adjacency(g), rank(g.n), label(g.n, 0), degeneracy(0),
      minCliqueSize(std::max<ui>(1, outputThreshold)), cliqueCount(0),
      maxCliqueSize(0), checksCount(0),
      hybridReorderSibling(useHybridReorderSibling),
      enableAdvancedRules(false), enableTailKernels(false),
      enableLocalBitset(false),
      siblingEventBudget(g.n < 50000 ? 8 : 64),
      siblingEvents(0), siblingBranchesBefore(0), siblingBranchesAfter(0),
      tinyKernelCalls(0), localBitsetHandoffs(0), localBitsetChecks(0),
      plex3Terminals(0), plex3Cliques(0), xDominanceRemoved(0),
      universalPForces(0), degreeZeroTerminals(0), degreeOneTerminals(0) {}

void FastListBK::emitClique(const std::vector<ui> &extension) const {
  if (!cliqueSink)
    return;

  std::vector<ui> clique = cliqueStack;
  clique.insert(clique.end(), extension.begin(), extension.end());
  std::sort(clique.begin(), clique.end());
  cliqueSink(clique);
}

void FastListBK::emitComplementMatching(const std::vector<ui> &p) const {
  if (!cliqueSink)
    return;

  std::vector<ui> forced;
  std::vector<std::pair<ui, ui>> missingEdges;
  std::vector<unsigned char> paired(p.size(), 0);
  for (size_t i = 0; i < p.size(); ++i) {
    if (paired[i] != 0)
      continue;
    size_t mate = p.size();
    for (size_t j = i + 1; j < p.size(); ++j) {
      if (paired[j] == 0 && !adjacency.contains(p[i], p[j])) {
        mate = j;
        break;
      }
    }
    if (mate == p.size()) {
      forced.push_back(p[i]);
    } else {
      paired[i] = paired[mate] = 1;
      missingEdges.push_back({p[i], p[mate]});
    }
  }

  std::vector<ui> extension = forced;
  std::function<void(size_t)> materialize = [&](size_t at) {
    if (at == missingEdges.size()) {
      emitClique(extension);
      return;
    }
    extension.push_back(missingEdges[at].first);
    materialize(at + 1);
    extension.back() = missingEdges[at].second;
    materialize(at + 1);
    extension.pop_back();
  };
  materialize(0);
}

void FastListBK::buildDegeneracyOrder() {
  std::vector<ui> degree(graph.degree.begin(), graph.degree.end());
  degeneracy = 0;
  ui maxDegree = 0;
  for (ui d : degree)
    maxDegree = std::max(maxDegree, d);

  if (graph.n != 0) {
    // Matula-Beck bin peeling: O(n + m), with no mutation of Graph::degree.
    std::vector<ui> bin(maxDegree + 1, 0);
    std::vector<ui> position(graph.n);
    std::vector<ui> vertices(graph.n);
    for (ui d : degree)
      ++bin[d];

    ui start = 0;
    for (ui d = 0; d <= maxDegree; ++d) {
      const ui count = bin[d];
      bin[d] = start;
      start += count;
    }
    for (ui u = 0; u < graph.n; ++u) {
      position[u] = bin[degree[u]]++;
      vertices[position[u]] = u;
    }
    for (ui d = maxDegree; d > 0; --d)
      bin[d] = bin[d - 1];
    bin[0] = 0;

    for (ui i = 0; i < graph.n; ++i) {
      const ui u = vertices[i];
      rank[u] = i;
      degeneracy = std::max(degeneracy, degree[u]);
      for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at) {
        const ui v = graph.neighbors[at];
        if (degree[v] > degree[u]) {
          const ui dv = degree[v];
          const ui pv = position[v];
          const ui pw = bin[dv];
          const ui w = vertices[pw];
          if (v != w) {
            position[v] = pw;
            position[w] = pv;
            vertices[pv] = w;
            vertices[pw] = v;
          }
          ++bin[dv];
          --degree[v];
        }
      }
    }
  }

  levels.clear();
  levels.resize(static_cast<size_t>(degeneracy) + 2);
#ifdef FASTLIST_OPPORTUNITY_PROFILE
  profilePAtDepth.assign(levels.size(), 0);
  profileStateAtDepth.assign(levels.size(), 0);
#endif
  if (levels.size() > 1) {
    levels[1].p.reserve(degeneracy);
    levels[1].x.reserve(maxDegree);
    levels[1].branch.reserve(degeneracy);
  }
}

ui FastListBK::neighborsInP(ui u, ui depth, const std::vector<ui> &p,
                            bool haveIncumbent, ui incumbent,
                            bool candidateFromX) {
  ui count = 0;
  const bool hashPath = graph.degree[u] > p.size();
  // Scores of P vertices feed the exact clique/plex/reduction tests below and
  // must remain exact. For an X candidate, however, once its best possible
  // score cannot exceed the incumbent it cannot change either the pivot or the
  // X-universal test, so the rest of the intersection is unnecessary.
  // Short scans stay on the branch-free baseline loop: checking a bound after
  // every item costs more than it saves there.
  constexpr ull EARLY_EXIT_MIN_ITEMS = 16;
  const ull potentialItems =
      hashPath ? static_cast<ull>(p.size()) : graph.degree[u];
  const bool allowEarlyExit =
#ifdef FASTLIST_DISABLE_PIVOT_EARLY_EXIT
      false;
#else
      enableAdvancedRules && candidateFromX && haveIncumbent &&
          potentialItems >= EARLY_EXIT_MIN_ITEMS;
#endif
#ifdef FASTLIST_OPPORTUNITY_PROFILE
  const ull scanItems = hashPath ? static_cast<ull>(p.size())
                                 : static_cast<ull>(graph.degree[u]);
  if (candidateFromX) {
    ++profile.pivotCandidatesX;
    if (hashPath)
      profile.pivotHashItemsX += scanItems;
    else
      profile.pivotCsrItemsX += scanItems;
  } else {
    ++profile.pivotCandidatesP;
    if (hashPath)
      profile.pivotHashItemsP += scanItems;
    else
      profile.pivotCsrItemsP += scanItems;
  }
  profile.pivotItems[profileBucket(profileStateAtDepth[depth])] += scanItems;

  bool abortable = false;
  ull savedItems = 0;
  if (allowEarlyExit && !hashPath &&
      std::min<ull>(graph.degree[u], p.size()) <= incumbent) {
    abortable = true;
    savedItems = scanItems;
  }
#endif

  if (allowEarlyExit && !hashPath &&
      std::min<ull>(graph.degree[u], p.size()) <= incumbent) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    ++profile.pivotAbortableX;
    profile.pivotAbortSavedX += scanItems;
#endif
    // Returning the incumbent preserves the strict score comparison while
    // avoiding a special "inexact" sentinel in the hot pivot loop.
    return incumbent;
  }

  if (hashPath) {
    if (allowEarlyExit) {
      size_t seen = 0;
      for (ui v : p) {
        count += adjacency.contains(u, v);
        ++seen;
        const ull remaining = p.size() - seen;
        if (static_cast<ull>(count) + remaining > incumbent)
          continue;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
        abortable = true;
        savedItems = remaining;
#endif
        break;
      }
    } else {
      for (ui v : p)
        count += adjacency.contains(u, v);
    }
  } else {
    const ui end = graph.offset[u + 1];
    if (allowEarlyExit) {
      for (ui at = graph.offset[u]; at < end; ++at) {
        count += label[graph.neighbors[at]] == static_cast<int>(depth);
        const ull remaining = end - at - 1;
        const ull possibleP =
            p.size() > count ? static_cast<ull>(p.size() - count) : 0;
        if (static_cast<ull>(count) + std::min(remaining, possibleP) >
            incumbent)
          continue;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
        abortable = true;
        savedItems = remaining;
#endif
        break;
      }
    } else {
      for (ui at = graph.offset[u]; at < end; ++at)
        count += label[graph.neighbors[at]] == static_cast<int>(depth);
    }
  }

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  if (abortable) {
    if (candidateFromX) {
      ++profile.pivotAbortableX;
      profile.pivotAbortSavedX += savedItems;
    } else {
      ++profile.pivotAbortableP;
      profile.pivotAbortSavedP += savedItems;
    }
  }
#endif
  return count;
}

ui FastListBK::neighborsInPBaseline(ui u, ui depth,
                                    const std::vector<ui> &p) const {
  ui count = 0;
  if (graph.degree[u] > p.size()) {
    for (ui v : p)
      count += adjacency.contains(u, v);
  } else {
    for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at)
      count += label[graph.neighbors[at]] == static_cast<int>(depth);
  }
  return count;
}

bool FastListBK::solveTinyP(ui cliqueSize, Level &level) {
  ++tinyKernelCalls;
  const ui pSize = static_cast<ui>(level.p.size());
  const ui subsetCount = 1U << pSize;
  std::array<ui, 1U << TINY_P_LIMIT> pAdjacency{};
  std::array<bool, 1U << TINY_P_LIMIT> blockedByX{};

  // Encode the induced graph on P. At this size, enumerating its subsets is
  // cheaper than constructing child vectors and revisiting the same edges.
  for (ui i = 0; i < pSize; ++i) {
    ui mask = 0;
    for (ui j = 0; j < pSize; ++j) {
      if (i != j && adjacency.contains(level.p[i], level.p[j]))
        mask |= 1U << j;
    }
    pAdjacency[i] = mask;
  }

  // Record the P-neighbourhood of every forbidden vertex, then use a small
  // superset zeta transform: blockedByX[Q] is true iff some x is adjacent to
  // every vertex of Q.
  for (ui x : level.x) {
    ui mask = 0;
    for (ui i = 0; i < pSize; ++i) {
      if (adjacency.contains(x, level.p[i]))
        mask |= 1U << i;
    }
    blockedByX[mask] = true;
  }
  for (ui bit = 0; bit < pSize; ++bit) {
    for (ui mask = 0; mask < subsetCount; ++mask) {
      if ((mask & (1U << bit)) == 0)
        blockedByX[mask] =
            blockedByX[mask] || blockedByX[mask | (1U << bit)];
    }
  }

  bool foundAny = false;
  ull foundHere = 0;
  for (ui subset = 1; subset < subsetCount; ++subset) {
    bool isClique = true;
    for (ui i = 0; i < pSize; ++i) {
      const ui bit = 1U << i;
      if ((subset & bit) != 0 &&
          (pAdjacency[i] & (subset ^ bit)) != (subset ^ bit)) {
        isClique = false;
        break;
      }
    }
    if (!isClique || blockedByX[subset])
      continue;

    bool maximalInP = true;
    for (ui i = 0; i < pSize; ++i) {
      if ((subset & (1U << i)) == 0 &&
          (pAdjacency[i] & subset) == subset) {
        maximalInP = false;
        break;
      }
    }
    if (!maximalInP)
      continue;

    const ui maximalSize =
        cliqueSize + static_cast<ui>(__builtin_popcount(subset));
    if (maximalSize < minCliqueSize)
      continue;

    addCliqueCountOrThrow(cliqueCount, 1);
    ++foundHere;
    maxCliqueSize = std::max(maxCliqueSize, maximalSize);
    if (cliqueSink) {
      std::vector<ui> extension;
      extension.reserve(__builtin_popcount(subset));
      for (ui i = 0; i < pSize; ++i)
        if ((subset & (1U << i)) != 0)
          extension.push_back(level.p[i]);
      emitClique(extension);
    }
    foundAny = true;
    if (hybridReorderSibling && siblingEvents < siblingEventBudget &&
        level.witness.empty()) {
      level.witness = cliqueStack;
      for (ui i = 0; i < pSize; ++i) {
        if ((subset & (1U << i)) != 0)
          level.witness.push_back(level.p[i]);
      }
    }
  }

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  ++profile.tinyKernelNodes;
  addCliqueCountOrThrow(profile.tinyKernelCliques, foundHere);
#endif
  return foundAny;
}

void FastListBK::reduceDominatedX(ui depth, Level &level) {
  const size_t pSize = level.p.size();
  const size_t xSize = level.x.size();
  // Dominance is quadratic in X and duplicates work done by the word kernel
  // on tiny states. Restrict it to a small forbidden frontier paired with a
  // substantial P, where removing an X vertex can save many descendant
  // maximality and pivot intersections.
  if (pSize < 16 || xSize < 2 || xSize > 8 || pSize + xSize > 64)
    return;

  std::vector<ull> masks(xSize, 0);
  for (size_t xi = 0; xi < xSize; ++xi) {
    for (size_t pi = 0; pi < pSize; ++pi) {
      if (adjacency.contains(level.x[xi], level.p[pi]))
        masks[xi] |= 1ULL << pi;
    }
  }

  std::vector<ui> active;
  active.reserve(xSize);
  const int restoredLabel = depth == 1 ? 0 : -static_cast<int>(depth - 1);
  for (size_t i = 0; i < xSize; ++i) {
    // An X vertex with no P neighbor cannot extend any continuation from this
    // nonempty-P state. For nonempty masks retain only inclusion-maximal
    // P-neighborhoods; equal masks keep their first representative.
    bool dominated = masks[i] == 0;
    const int popI = __builtin_popcountll(masks[i]);
    for (size_t j = 0; !dominated && j < xSize; ++j) {
      if (i == j || __builtin_popcountll(masks[j]) < popI)
        continue;
      dominated = (masks[i] & ~masks[j]) == 0 &&
                  (masks[i] != masks[j] || j < i);
    }
    if (dominated) {
      label[level.x[i]] = restoredLabel;
      ++xDominanceRemoved;
    } else {
      active.push_back(level.x[i]);
    }
  }
  level.x.swap(active);
}

bool FastListBK::solveLowDegreeChild(ui cliqueSize, Level &child,
                                     bool &found) {
  if (child.p.size() > 1)
    return false;

  incrementSearchStateOrThrow(checksCount);
  child.witness.clear();
  found = false;
  if (child.p.empty()) {
    ++degreeZeroTerminals;
    if (!child.x.empty() || cliqueSize < minCliqueSize)
      return true;
    addCliqueCountOrThrow(cliqueCount, 1);
    maxCliqueSize = std::max(maxCliqueSize, cliqueSize);
    if (cliqueSink)
      emitClique();
    found = true;
    if (hybridReorderSibling && siblingEvents < siblingEventBudget)
      child.witness = cliqueStack;
    return true;
  }

  ++degreeOneTerminals;
  const ui extension = child.p.front();
  for (ui x : child.x) {
    if (adjacency.contains(x, extension))
      return true;
  }

  const ui maximalSize = cliqueSize + 1;
  if (maximalSize < minCliqueSize)
    return true;
  addCliqueCountOrThrow(cliqueCount, 1);
  maxCliqueSize = std::max(maxCliqueSize, maximalSize);
  if (cliqueSink)
    emitClique({extension});
  found = true;
  if (hybridReorderSibling && siblingEvents < siblingEventBudget) {
    child.witness = cliqueStack;
    child.witness.push_back(extension);
  }
  return true;
}

void FastListBK::intersectInto(ui u, ui depth, const Level &parent,
                               Level &child) {
  child.p.clear();
  child.x.clear();
  const int pLabel = static_cast<int>(depth);
  const int childLabel = static_cast<int>(depth + 1);

  if (graph.degree[u] > parent.p.size() + parent.x.size()) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    const ull items = parent.p.size() + parent.x.size();
    ++profile.intersectHashCalls;
    profile.intersectHashItems += items;
    profile.intersectItems[profileBucket(profileStateAtDepth[depth])] += items;
#endif
    for (ui v : parent.p) {
      if (!adjacency.contains(u, v))
        continue;
      if (label[v] == pLabel) {
        child.p.push_back(v);
        label[v] = childLabel;
      } else if (label[v] == -pLabel) {
        child.x.push_back(v);
        label[v] = -childLabel;
      }
    }
    for (ui v : parent.x) {
      if (label[v] == -pLabel && adjacency.contains(u, v)) {
        child.x.push_back(v);
        label[v] = -childLabel;
      }
    }
  } else {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    const ull items = graph.degree[u];
    ++profile.intersectCsrCalls;
    profile.intersectCsrItems += items;
    profile.intersectItems[profileBucket(profileStateAtDepth[depth])] += items;
#endif
    for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at) {
      const ui v = graph.neighbors[at];
      if (label[v] == pLabel) {
        child.p.push_back(v);
        label[v] = childLabel;
      } else if (label[v] == -pLabel) {
        child.x.push_back(v);
        label[v] = -childLabel;
      }
    }
  }
}

bool FastListBK::trySiblingEffect(ui depth, size_t nextBranch,
                                  const std::vector<ui> &witness) {
  if (!hybridReorderSibling || siblingEvents >= siblingEventBudget ||
      witness.empty())
    return false;

  Level &level = levels[depth];
  const size_t remaining = level.branch.size() - nextBranch;
  if (remaining == 0)
    return false;

  // C is an emitted maximal clique containing the current R. Every distinct
  // maximal continuation must therefore contain a live vertex outside C.
  // Branching sequentially on liveP \\ C is the one-covering-clique sibling
  // effect; C-intersection vertices stay in P and remain available to children.
  std::vector<ui> sortedWitness = witness;
  std::sort(sortedWitness.begin(), sortedWitness.end());
  std::vector<ui> siblingRoots;
  siblingRoots.reserve(remaining);
  const int liveLabel = static_cast<int>(depth);
  for (ui v : level.p) {
    if (label[v] == liveLabel &&
        !std::binary_search(sortedWitness.begin(), sortedWitness.end(), v))
      siblingRoots.push_back(v);
  }

  // Use the sibling effect only when it does not enlarge the current pivot
  // frontier. Other nodes remain ordinary pivot subbranches.
  if (siblingRoots.size() > remaining)
    return false;
  if (siblingRoots.size() == remaining &&
      std::equal(siblingRoots.begin(), siblingRoots.end(),
                 level.branch.begin() + nextBranch))
    return false;

  siblingBranchesBefore += remaining;
  siblingBranchesAfter += siblingRoots.size();
  ++siblingEvents;
  level.branch.swap(siblingRoots);
  return true;
}

bool FastListBK::enumerateBaseline(ui depth, ui cliqueSize) {
  incrementSearchStateOrThrow(checksCount);
  Level &level = levels[depth];
  level.witness.clear();
  if (level.p.empty()) {
    if (level.x.empty() && cliqueSize >= minCliqueSize) {
      addCliqueCountOrThrow(cliqueCount, 1);
      maxCliqueSize = std::max(maxCliqueSize, cliqueSize);
      if (cliqueSink)
        emitClique();
      if (hybridReorderSibling && siblingEvents < siblingEventBudget)
        level.witness = cliqueStack;
      return true;
    }
    return false;
  }

  ui pivot = level.x.empty() ? level.p.front() : level.x.front();
  ui best = 0;
  bool havePivot = false;
  const ui pSize = static_cast<ui>(level.p.size());
  for (ui u : level.x) {
    const ui score = neighborsInPBaseline(u, depth, level.p);
    // An X vertex adjacent to all of P proves every continuation non-maximal.
    if (score == pSize)
      return false;
    if (!havePivot || score > best) {
      pivot = u;
      best = score;
      havePivot = true;
    }
  }

  ui minPScore = pSize;
  ui deficientByOne = 0;
  for (ui u : level.p) {
    const ui score = neighborsInPBaseline(u, depth, level.p);
    minPScore = std::min(minPScore, score);
    if (pSize >= 2 && score + 2 == pSize)
      ++deficientByOne;
    if (!havePivot || score > best) {
      pivot = u;
      best = score;
      havePivot = true;
    }
  }

  // P is a clique. Since no X vertex was universal, R union P is the one
  // maximal continuation from this state.
  if (minPScore + 1 == pSize) {
    const ui maximalSize = cliqueSize + pSize;
    if (maximalSize >= minCliqueSize) {
      addCliqueCountOrThrow(cliqueCount, 1);
      maxCliqueSize = std::max(maxCliqueSize, maximalSize);
      if (cliqueSink)
        emitClique(level.p);
      if (hybridReorderSibling && siblingEvents < siblingEventBudget) {
        level.witness = cliqueStack;
        level.witness.insert(level.witness.end(), level.p.begin(),
                             level.p.end());
      }
      return true;
    }
    return false;
  }

  // With X empty and complement degree at most one, the complement of P is a
  // matching. Each missing-edge pair contributes an independent endpoint
  // choice, so there are exactly 2^k maximal continuations of equal size.
  if (level.x.empty() && pSize >= 2 && minPScore + 2 >= pSize &&
      (deficientByOne & 1U) == 0) {
    const ui missingEdges = deficientByOne >> 1;
    const ui maximalSize = cliqueSize + pSize - missingEdges;
    if (maximalSize < minCliqueSize)
      return false;
    if (missingEdges >= 64)
      throw std::overflow_error(
          "complement-matching clique count exceeds uint64_t");
    const ull add = 1ULL << missingEdges;
    addCliqueCountOrThrow(cliqueCount, add);
    maxCliqueSize = std::max(maxCliqueSize, maximalSize);
    if (cliqueSink)
      emitComplementMatching(level.p);
    // The closed form counts several cliques. Do not drive sibling
    // reordering without materializing one proven representative.
    return true;
  }

  level.branch.clear();
  level.processedRoots.clear();
  for (ui u : level.p) {
    if (!adjacency.contains(pivot, u))
      level.branch.push_back(u);
  }

  bool foundAny = false;
  size_t nextBranch = 0;
  while (nextBranch < level.branch.size()) {
    const ui u = level.branch[nextBranch++];
    if (label[u] != static_cast<int>(depth))
      continue;
    Level &child = levels[depth + 1];
    intersectInto(u, depth, level, child);
    if (needsCliqueStack())
      cliqueStack.push_back(u);
    const bool childFound = enumerateBaseline(depth + 1, cliqueSize + 1);
    if (needsCliqueStack())
      cliqueStack.pop_back();

    const int parentLabel = static_cast<int>(depth);
    for (ui v : child.p)
      label[v] = parentLabel;
    for (ui v : child.x)
      label[v] = -parentLabel;
    label[u] = -parentLabel;
    level.processedRoots.push_back(u);

    if (childFound) {
      foundAny = true;
      if (level.witness.empty() && !child.witness.empty())
        level.witness = child.witness;

      if (!child.witness.empty() &&
          trySiblingEffect(depth, nextBranch, child.witness))
        nextBranch = 0;
    }
  }

  const int parentLabel = static_cast<int>(depth);
  for (ui u : level.processedRoots)
    label[u] = parentLabel;
  return foundAny;
}

bool FastListBK::enumerate(ui depth, ui cliqueSize) {
#ifndef FASTLIST_OPPORTUNITY_PROFILE
  // In the local-only portfolio, |P| can only shrink below this point. Once
  // it is too small for the adaptive word kernel, keep the whole remaining
  // subtree on the zero-overhead list recursion as well.
#if defined(FASTLIST_DISABLE_LOCAL_BITSET) ||                              \
    defined(FASTLIST_DISABLE_LOCAL_ADAPTIVE)
  constexpr bool localAdaptiveCompiled = false;
#else
  constexpr bool localAdaptiveCompiled = true;
#endif
  const bool canReachLocalBitset =
      localAdaptiveCompiled && enableLocalBitset &&
      levels[depth].p.size() >= 12;
  if (!enableAdvancedRules && !enableTailKernels && !canReachLocalBitset)
    return enumerateBaseline(depth, cliqueSize);
#endif
  incrementSearchStateOrThrow(checksCount);
  Level &level = levels[depth];
  level.witness.clear();

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  const ui entryP = static_cast<ui>(level.p.size());
  const ui entryX = static_cast<ui>(level.x.size());
  const ui entryState = entryP + entryX;
  const ui parentP =
      depth == 1 ? std::numeric_limits<ui>::max() : profilePAtDepth[depth - 1];
  const ui parentState = depth == 1 ? std::numeric_limits<ui>::max()
                                    : profileStateAtDepth[depth - 1];
  profilePAtDepth[depth] = entryP;
  profileStateAtDepth[depth] = entryState;
  ++profile.stateNodes[profileBucket(entryState)];
  const size_t tinyP = std::min<size_t>(entryP, 16);
  const size_t tinyX = std::min<size_t>(entryX, 16);
  ++profile.tinyPX[tinyP * 17 + tinyX];

  for (size_t i = 0; i < PROFILE_THRESHOLD_COUNT; ++i) {
    const ui threshold = PROFILE_THRESHOLDS[i];
    if (entryState <= threshold) {
      ++profile.fullEligible[i];
      if (parentState > threshold) {
        if (entryP == 0) {
          ++profile.fullEmptyCrossings[i];
        } else {
          ++profile.fullFrontiers[i];
          profile.fullFrontierP[i] += entryP;
          profile.fullFrontierX[i] += entryX;
          profile.fullCells[i] += static_cast<ull>(entryState) * entryP;
          profile.fullWords[i] +=
              static_cast<ull>(entryState) * ((entryP + 63) >> 6);
        }
      }
    }
    if (entryP <= threshold) {
      ++profile.partialEligible[i];
      if (parentP > threshold) {
        if (entryP == 0) {
          ++profile.partialEmptyCrossings[i];
        } else {
          ++profile.partialFrontiers[i];
          profile.partialFrontierP[i] += entryP;
          profile.partialFrontierX[i] += entryX;
          profile.partialCells[i] += static_cast<ull>(entryP) * entryP;
          profile.partialWords[i] +=
              static_cast<ull>(entryP) * ((entryP + 63) >> 6);
        }
      }
    }
  }

  bool profilePlex3Started = false;
  ull profilePlex3ChecksStart = 0;
  ull profilePlex3CliquesStart = 0;
  ull profilePlex3PivotStart = 0;
  ull profilePlex3IntersectStart = 0;
  auto profileFinish = [&](bool result) {
    if (profilePlex3Started) {
      profile.plex3DescendantChecks +=
          checksCount - profilePlex3ChecksStart;
      addCliqueCountOrThrow(profile.plex3DescendantCliques,
                            cliqueCount - profilePlex3CliquesStart);
      profile.plex3DescendantPivotItems +=
          profilePivotItemsTotal() - profilePlex3PivotStart;
      profile.plex3SubtreeIntersectItems +=
          profileIntersectItemsTotal() - profilePlex3IntersectStart;
      profilePlex3Active = false;
    }
    return result;
  };
#define FASTLIST_PROFILE_RETURN(value) return profileFinish(value)
#else
#define FASTLIST_PROFILE_RETURN(value) return value
#endif

  if (level.p.empty()) {
    if (level.x.empty() && cliqueSize >= minCliqueSize) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      ++profile.leafMaximal;
#endif
      addCliqueCountOrThrow(cliqueCount, 1);
      maxCliqueSize = std::max(maxCliqueSize, cliqueSize);
      if (cliqueSink)
        emitClique();
      if (hybridReorderSibling && siblingEvents < siblingEventBudget)
        level.witness = cliqueStack;
      FASTLIST_PROFILE_RETURN(true);
    }
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    ++profile.leafBlockedX;
#endif
    FASTLIST_PROFILE_RETURN(false);
  }

  // A complete exact local kernel avoids pivot selection, child-vector
  // construction, and recursion at the overwhelmingly common tail states.
#ifndef FASTLIST_DISABLE_TINY_KERNEL
  if (enableTailKernels && level.p.size() <= TINY_P_LIMIT &&
      level.p.size() + level.x.size() <= TINY_P_LIMIT)
    FASTLIST_PROFILE_RETURN(solveTinyP(cliqueSize, level));
#endif

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  // Preserve the profiler's view of the unreduced state.
  profileXDominance(depth, level);
#endif
#ifndef FASTLIST_DISABLE_X_DOMINANCE
  if (enableAdvancedRules)
    reduceDominatedX(depth, level);
#endif

  ui pivot = level.x.empty() ? level.p.front() : level.x.front();
  ui best = 0;
  bool havePivot = false;
  const ui pSize = static_cast<ui>(level.p.size());
  ui universalCandidate = std::numeric_limits<ui>::max();
#ifdef FASTLIST_OPPORTUNITY_PROFILE
  bool pivotFromX = false;
  bool naudeWidth2Seen = false;
  ull naudeTailCandidates = 0;
  ull naudeTailItems = 0;
  ui pDeg0 = 0;
  ui pDeg1 = 0;
  ui pUniversal = 0;
  ull complementDegreeSum = 0;
#endif
  for (ui u : level.x) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    const ull candidateItemsBefore = profilePivotItemsTotal();
#endif
    const ui score =
        neighborsInP(u, depth, level.p, havePivot, best, true);
    // An X vertex adjacent to all of P proves every continuation non-maximal.
    if (score == pSize) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      ++profile.xUniversalPrune;
#endif
      FASTLIST_PROFILE_RETURN(false);
    }
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    const ull candidateItems =
        profilePivotItemsTotal() - candidateItemsBefore;
    if (naudeWidth2Seen) {
      ++naudeTailCandidates;
      naudeTailItems += candidateItems;
    } else if (pSize - score <= 2) {
      naudeWidth2Seen = true;
    }
#endif
    if (!havePivot || score > best) {
      pivot = u;
      best = score;
      havePivot = true;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      pivotFromX = true;
#endif
    }
  }

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  const bool hadXPivot = havePivot;
  const ui bestAfterX = best;
#endif
  ui minPScore = pSize;
  ui deficientByOne = 0;
  for (ui u : level.p) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    const ull candidateItemsBefore = profilePivotItemsTotal();
#endif
    const ui score =
        neighborsInP(u, depth, level.p, havePivot, best, false);
    minPScore = std::min(minPScore, score);
    if (pSize >= 2 && score + 2 == pSize)
      ++deficientByOne;
    if (enableAdvancedRules && score + 1 == pSize &&
        universalCandidate == std::numeric_limits<ui>::max())
      universalCandidate = u;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    pDeg0 += score == 0;
    pDeg1 += score == 1;
    pUniversal += score + 1 == pSize;
    complementDegreeSum += pSize - 1 - score;
    const ull candidateItems =
        profilePivotItemsTotal() - candidateItemsBefore;
    if (naudeWidth2Seen) {
      ++naudeTailCandidates;
      naudeTailItems += candidateItems;
    } else if (pSize - score <= 2) {
      naudeWidth2Seen = true;
    }
#endif
    if (!havePivot || score > best) {
      pivot = u;
      best = score;
      havePivot = true;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      pivotFromX = false;
#endif
    }
  }

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  if (naudeWidth2Seen) {
    ++profile.naudeWidth2Nodes;
    profile.naudeTailCandidates += naudeTailCandidates;
    profile.naudeTailItems += naudeTailItems;
  }
  profile.pImprovedPivotNodes += hadXPivot && best > bestAfterX;
  profile.pivotWinnerX += pivotFromX;
  profile.pivotWinnerP += !pivotFromX;
#endif

  // P is a clique. Since no X vertex was universal, R union P is the one
  // maximal continuation from this state.
  if (minPScore + 1 == pSize) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    ++profile.pCliqueSolved;
#endif
    const ui maximalSize = cliqueSize + pSize;
    if (maximalSize >= minCliqueSize) {
      addCliqueCountOrThrow(cliqueCount, 1);
      maxCliqueSize = std::max(maxCliqueSize, maximalSize);
      if (cliqueSink)
        emitClique(level.p);
      if (hybridReorderSibling && siblingEvents < siblingEventBudget) {
        level.witness = cliqueStack;
        level.witness.insert(level.witness.end(), level.p.begin(),
                             level.p.end());
      }
      FASTLIST_PROFILE_RETURN(true);
    }
    FASTLIST_PROFILE_RETURN(false);
  }

  // With X empty and complement degree at most one, the complement of P is a
  // matching. Each missing-edge pair contributes an independent endpoint
  // choice, so there are exactly 2^k maximal continuations of equal size.
  if (level.x.empty() && pSize >= 2 && minPScore + 2 >= pSize &&
      (deficientByOne & 1U) == 0) {
    const ui missingEdges = deficientByOne >> 1;
    const ui maximalSize = cliqueSize + pSize - missingEdges;
    if (maximalSize < minCliqueSize) {
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      ++profile.matchingSolved;
#endif
      FASTLIST_PROFILE_RETURN(false);
    }
    if (missingEdges < 64) {
      const ull add = 1ULL << missingEdges;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
      ++profile.matchingSolved;
#endif
      addCliqueCountOrThrow(cliqueCount, add);
      maxCliqueSize = std::max(maxCliqueSize, maximalSize);
      if (cliqueSink)
        emitComplementMatching(level.p);
      // The closed form counts several cliques. Do not drive sibling
      // reordering without materializing one proven representative.
      FASTLIST_PROFILE_RETURN(true);
    }
#ifdef FASTLIST_OPPORTUNITY_PROFILE
    ++profile.matchingOverflowFallback;
#endif
    throw std::overflow_error(
        "complement-matching clique count exceeds uint64_t");
  }

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  if (pDeg0 != 0) {
    ++profile.pDeg0Nodes;
    profile.pDeg0Vertices += pDeg0;
    if (entryState > 3)
      profile.pDeg0OutsideTiny3 += pDeg0;
  }
  if (pDeg1 != 0) {
    ++profile.pDeg1Nodes;
    profile.pDeg1Vertices += pDeg1;
    if (entryState > 3)
      profile.pDeg1OutsideTiny3 += pDeg1;
  }
  if (pUniversal != 0) {
    ++profile.universalPNodes;
    profile.universalPVertices += pUniversal;
    if (entryState > 3)
      profile.universalPOutsideTiny3 += pUniversal;
  }

  if (level.x.empty() && minPScore + 3 == pSize) {
    ++profile.plex3Nodes;
    profile.plex3PVertices += pSize;
    profile.plex3ComplementEdges += complementDegreeSum >> 1;
    if (!profilePlex3Active) {
      profilePlex3Active = true;
      profilePlex3Started = true;
      ++profile.plex3TopmostNodes;
      profilePlex3ChecksStart = checksCount;
      profilePlex3CliquesStart = cliqueCount;
      profilePlex3PivotStart = profilePivotItemsTotal();
      profilePlex3IntersectStart = profileIntersectItemsTotal();
    }
  }

#endif

  // HBBMC's 3-plex terminal: with X empty and complement degree at most two,
  // maximal clique continuations are maximal independent sets of disjoint
  // complement paths/cycles and can be counted without BK branching.
#ifndef FASTLIST_DISABLE_PLEX3
  if (enableAdvancedRules && level.x.empty() && minPScore + 3 >= pSize) {
    std::vector<std::vector<ui>> materializedCliques;
    FastCliqueSink bufferedSink;
    if (cliqueSink) {
      bufferedSink = [&](const std::vector<ui> &clique) {
        materializedCliques.push_back(clique);
      };
    }
    const std::vector<ui> *prefix =
        cliqueSink ||
                (hybridReorderSibling && siblingEvents < siblingEventBudget)
            ? &cliqueStack
            : nullptr;
    FastPlex3Result plex =
        solveFastPlex3Subtree(adjacency, level.p, cliqueSize, prefix,
                              minCliqueSize,
                              cliqueSink ? &bufferedSink : nullptr);
    if (plex.handled) {
      ull combinedCliqueCount = 0;
      ull combinedPlex3Cliques = 0;
      if (!tryAddUll(cliqueCount, plex.cliqueCount, combinedCliqueCount) ||
          !tryAddUll(plex3Cliques, plex.cliqueCount,
                     combinedPlex3Cliques))
        throw std::overflow_error(
            "3-plex aggregate clique count exceeds uint64_t");
      ++plex3Terminals;
      plex3Cliques = combinedPlex3Cliques;
      cliqueCount = combinedCliqueCount;
      maxCliqueSize = std::max(maxCliqueSize, plex.maxCliqueSize);
      if (cliqueSink)
        for (const std::vector<ui> &clique : materializedCliques)
          cliqueSink(clique);
      level.witness = std::move(plex.witness);
      FASTLIST_PROFILE_RETURN(plex.found);
    }
  }
#endif

  // A P-universal candidate belongs to every maximal continuation. Force it
  // into R and recurse once instead of generating sibling branches.
#ifndef FASTLIST_DISABLE_UNIVERSAL_P
  if (enableAdvancedRules &&
      universalCandidate != std::numeric_limits<ui>::max()) {
    ++universalPForces;
    Level &child = levels[depth + 1];
    intersectInto(universalCandidate, depth, level, child);
    if (needsCliqueStack())
      cliqueStack.push_back(universalCandidate);
    const bool childFound = enumerate(depth + 1, cliqueSize + 1);
    if (needsCliqueStack())
      cliqueStack.pop_back();

    const int parentLabel = static_cast<int>(depth);
    for (ui v : child.p)
      label[v] = parentLabel;
    for (ui v : child.x)
      label[v] = -parentLabel;
    if (childFound)
      level.witness = child.witness;
    FASTLIST_PROFILE_RETURN(childFound);
  }
#endif

  // Convert a small frontier once and consume its entire subtree with the
  // word-parallel kernel. The branch-width gate avoids matrix setup on easy
  // list states; very small full states use this as a second tiny kernel.
  const size_t stateSize = level.p.size() + level.x.size();
  const ui branchWidth = pSize - best;
  const bool tinyBitsetState =
#ifdef FASTLIST_DISABLE_LOCAL_TINY
      false;
#else
      enableTailKernels && stateSize <= 12;
#endif
  const bool adaptiveBitsetState =
#ifdef FASTLIST_DISABLE_LOCAL_ADAPTIVE
      false;
#else
      enableLocalBitset && pSize >= 12 && stateSize <= 64 &&
          branchWidth >= 4;
#endif
#ifndef FASTLIST_DISABLE_LOCAL_BITSET
  if (tinyBitsetState || adaptiveBitsetState) {
    std::vector<std::vector<ui>> materializedCliques;
    FastCliqueSink bufferedSink;
    if (cliqueSink) {
      bufferedSink = [&](const std::vector<ui> &clique) {
        materializedCliques.push_back(clique);
      };
    }
    const std::vector<ui> *prefix =
        cliqueSink ||
                (hybridReorderSibling && siblingEvents < siblingEventBudget)
            ? &cliqueStack
            : nullptr;
    FastLocalBitsetResult local = solveFastLocalBitsetSubtree(
        adjacency, level.p, level.x, cliqueSize, prefix, minCliqueSize,
        cliqueSink ? &bufferedSink : nullptr);
    const ull extraChecks =
        local.checksCount == 0 ? 0 : local.checksCount - 1;
    if (local.handled) {
      ull combinedCliqueCount = 0;
      if (!tryAddUll(cliqueCount, local.cliqueCount, combinedCliqueCount))
        throw std::overflow_error(
            "local-bitset aggregate clique count exceeds uint64_t");
      if (extraChecks > std::numeric_limits<ull>::max() - checksCount)
        throw std::overflow_error("search check count exceeds uint64_t");
      ++localBitsetHandoffs;
      if (local.checksCount >
          std::numeric_limits<ull>::max() - localBitsetChecks)
        throw std::overflow_error(
            "local-bitset check count exceeds uint64_t");
      localBitsetChecks += local.checksCount;
      checksCount += extraChecks;
      cliqueCount = combinedCliqueCount;
      maxCliqueSize = std::max(maxCliqueSize, local.maxCliqueSize);
      if (cliqueSink)
        for (const std::vector<ui> &clique : materializedCliques)
          cliqueSink(clique);
      level.witness = std::move(local.witness);
      FASTLIST_PROFILE_RETURN(local.found);
    }
  }
#endif

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  ++profile.ordinaryBranchNodes;
  profile.ordinaryTiny3Nodes += entryState <= 3;
#endif

  level.branch.clear();
  level.processedRoots.clear();
  for (ui u : level.p) {
    if (!adjacency.contains(pivot, u))
      level.branch.push_back(u);
  }
#ifdef FASTLIST_OPPORTUNITY_PROFILE
  profile.ordinaryBranchVertices += level.branch.size();
#endif

  bool foundAny = false;
  size_t nextBranch = 0;
  while (nextBranch < level.branch.size()) {
    const ui u = level.branch[nextBranch++];
    if (label[u] != static_cast<int>(depth))
      continue;
    Level &child = levels[depth + 1];
    intersectInto(u, depth, level, child);
    if (needsCliqueStack())
      cliqueStack.push_back(u);
    bool childFound = false;
#ifndef FASTLIST_DISABLE_LOW_DEGREE
    if (!enableTailKernels ||
        !solveLowDegreeChild(cliqueSize + 1, child, childFound))
      childFound = enumerate(depth + 1, cliqueSize + 1);
#else
    childFound = enumerate(depth + 1, cliqueSize + 1);
#endif
    if (needsCliqueStack())
      cliqueStack.pop_back();

    const int parentLabel = static_cast<int>(depth);
    for (ui v : child.p)
      label[v] = parentLabel;
    for (ui v : child.x)
      label[v] = -parentLabel;
    label[u] = -parentLabel;
    level.processedRoots.push_back(u);

    if (childFound) {
      foundAny = true;
      if (level.witness.empty() && !child.witness.empty())
        level.witness = child.witness;

      if (!child.witness.empty() &&
          trySiblingEffect(depth, nextBranch, child.witness))
        nextBranch = 0;
    }
  }

  const int parentLabel = static_cast<int>(depth);
  for (ui u : level.processedRoots)
    label[u] = parentLabel;
  FASTLIST_PROFILE_RETURN(foundAny);
#undef FASTLIST_PROFILE_RETURN
}

void FastListBK::findAllMaximalCliques(const std::string &outputLabel) {
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  siblingEvents = 0;
  siblingBranchesBefore = 0;
  siblingBranchesAfter = 0;
  tinyKernelCalls = 0;
  localBitsetHandoffs = 0;
  localBitsetChecks = 0;
  plex3Terminals = 0;
  plex3Cliques = 0;
  xDominanceRemoved = 0;
  universalPForces = 0;
  degreeZeroTerminals = 0;
  degreeOneTerminals = 0;
  cliqueStack.clear();
  std::fill(label.begin(), label.end(), 0);

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  profile = OpportunityProfile{};
  profilePlex3Active = false;
  profileGraphRules();
#endif

  const auto start = std::chrono::high_resolution_clock::now();
  buildDegeneracyOrder();

  const ull possibleEdges =
      graph.n < 2 ? 0 : static_cast<ull>(graph.n) * (graph.n - 1) / 2;
  const double averageDegree =
      graph.n == 0 ? 0.0 : (2.0 * static_cast<double>(graph.m)) / graph.n;
  const bool hardDenseGraph =
      graph.n <= 4096 && possibleEdges != 0 &&
      static_cast<ull>(graph.m) * 5 >= possibleEdges;
  const bool giantLowDegeneracyGraph =
      graph.n >= 1000000 && degeneracy <= 4;
  const bool sparseDenseCoreGraph =
      graph.n >= 10000 && graph.n <= 200000 && averageDegree >= 4.0 &&
      averageDegree <= 16.0 && degeneracy >= 24 &&
      static_cast<double>(degeneracy) >= averageDegree * 3.5;
  enableAdvancedRules = hardDenseGraph;
  enableTailKernels = hardDenseGraph || giantLowDegeneracyGraph;
  enableLocalBitset = enableTailKernels || sparseDenseCoreGraph;

  if (graph.n != 0) {
    Level &root = levels[1];
    for (ui u = 0; u < graph.n; ++u) {
      root.p.clear();
      root.x.clear();
      for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at) {
        const ui v = graph.neighbors[at];
        if (rank[v] < rank[u]) {
          root.x.push_back(v);
          label[v] = -1;
        } else {
          root.p.push_back(v);
          label[v] = 1;
        }
      }

      if (needsCliqueStack())
        cliqueStack.push_back(u);
#ifndef FASTLIST_OPPORTUNITY_PROFILE
#if defined(FASTLIST_DISABLE_LOCAL_BITSET) ||                              \
    defined(FASTLIST_DISABLE_LOCAL_ADAPTIVE)
      constexpr bool localAdaptiveCompiled = false;
#else
      constexpr bool localAdaptiveCompiled = true;
#endif
      const bool rootCanReachLocalBitset =
          localAdaptiveCompiled && enableLocalBitset && root.p.size() >= 12;
      if (!enableAdvancedRules && !enableTailKernels &&
          !rootCanReachLocalBitset)
        enumerateBaseline(1, 1);
      else
        enumerate(1, 1);
#else
      enumerate(1, 1);
#endif
      if (needsCliqueStack())
        cliqueStack.pop_back();
      for (ui at = graph.offset[u]; at < graph.offset[u + 1]; ++at)
        label[graph.neighbors[at]] = 0;
    }
  }

  const auto finish = std::chrono::high_resolution_clock::now();
  const double ms =
      std::chrono::duration<double, std::milli>(finish - start).count();
  std::cout << outputLabel << ": cliques=" << cliqueCount
            << "  maxSize=" << maxCliqueSize
            << "  minSize=" << minCliqueSize << "  checks=" << checksCount
            << "  degeneracy=" << degeneracy
            << "  portfolio="
            << (enableAdvancedRules
                    ? "dense"
                    : (enableTailKernels
                           ? "lowdeg"
                           : (enableLocalBitset ? "local" : "baseline")))
            << "  siblingEvents=" << siblingEvents
            << "  siblingBranches=" << siblingBranchesBefore << "->"
            << siblingBranchesAfter
            << "  tiny=" << tinyKernelCalls
            << "  local=" << localBitsetHandoffs
            << "  plex3=" << plex3Terminals
            << "  xdom=" << xDominanceRemoved
            << "  universal=" << universalPForces
            << "  lowDegree=" << degreeZeroTerminals << "/"
            << degreeOneTerminals
            << "  time=" << std::fixed << std::setprecision(3) << ms << " ms"
            << std::endl;
#ifdef FASTLIST_OPPORTUNITY_PROFILE
  printOpportunityProfile();
#endif
}
