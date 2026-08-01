#include "fast_local_bitset.h"
#include "checked_count.h"

#include <array>
#include <limits>

namespace {

class OneWordSubtreeSolver {
private:
  const std::vector<ui> &vertices;
  const std::vector<ui> *cliquePrefix;
  const FastCliqueSink *cliqueSink;
  ui minCliqueSize;
  std::array<ull, 64> neighbors{};
  FastLocalBitsetResult result;
  std::vector<std::vector<ui>> pendingCliques;
  bool overflow = false;

  static ui popcount(ull bits) {
    return static_cast<ui>(__builtin_popcountll(bits));
  }

  bool addCheck() {
    if (result.checksCount == std::numeric_limits<ull>::max()) {
      overflow = true;
      return false;
    }
    ++result.checksCount;
    return true;
  }

  bool cliqueSizeWithExtension(ui cliqueSize, ui extensionSize,
                               ui &total) {
    if (extensionSize > std::numeric_limits<ui>::max() - cliqueSize) {
      overflow = true;
      return false;
    }
    total = cliqueSize + extensionSize;
    return true;
  }

  bool addCliques(ull count, ui size) {
    ull combined = 0;
    if (!tryAddUll(result.cliqueCount, count, combined)) {
      overflow = true;
      return false;
    }
    result.cliqueCount = combined;
    result.maxCliqueSize = std::max(result.maxCliqueSize, size);
    result.found = true;
    return true;
  }

  void saveWitness(ull selected) {
    if (cliquePrefix == nullptr || !result.witness.empty())
      return;

    result.witness = *cliquePrefix;
    while (selected != 0) {
      const ui index = static_cast<ui>(__builtin_ctzll(selected));
      result.witness.push_back(vertices[index]);
      selected &= selected - 1;
    }
  }

  void materialize(ull selected) {
    if (cliqueSink == nullptr)
      return;

    std::vector<ui> clique;
    if (cliquePrefix != nullptr)
      clique = *cliquePrefix;
    while (selected != 0) {
      const ui index = static_cast<ui>(__builtin_ctzll(selected));
      clique.push_back(vertices[index]);
      selected &= selected - 1;
    }
    std::sort(clique.begin(), clique.end());
    pendingCliques.push_back(std::move(clique));
  }

  void materializeMatching(ull selected, ull remaining) {
    if (remaining == 0) {
      materialize(selected);
      return;
    }

    const ui index = static_cast<ui>(__builtin_ctzll(remaining));
    const ull bit = 1ULL << index;
    remaining &= ~bit;
    const ull nonNeighbors = remaining & ~neighbors[index];
    if (nonNeighbors == 0) {
      materializeMatching(selected | bit, remaining);
      return;
    }

    const ull mate = nonNeighbors & (~nonNeighbors + 1);
    remaining &= ~mate;
    materializeMatching(selected | bit, remaining);
    materializeMatching(selected | mate, remaining);
  }

  // When the complement of P is a matching, every maximal clique contains
  // every complement-isolated vertex and one endpoint of every missing edge.
  // Return one such clique for sibling reordering while the closed form counts
  // all of them.
  ull matchingWitness(ull p) const {
    ull selected = 0;
    ull remaining = p;
    while (remaining != 0) {
      const ui index = static_cast<ui>(__builtin_ctzll(remaining));
      const ull bit = 1ULL << index;
      remaining &= ~bit;

      const ull nonNeighbors = remaining & ~neighbors[index];
      selected |= bit;
      if (nonNeighbors != 0) {
        // The caller has already established complement degree at most one.
        remaining &= ~nonNeighbors;
      }
    }
    return selected;
  }

  void enumerate(ull p, ull x, ui cliqueSize, ull chosen) {
    if (overflow || !addCheck())
      return;

    if (p == 0) {
      if (x == 0 && cliqueSize >= minCliqueSize) {
        if (addCliques(1, cliqueSize)) {
          saveWitness(chosen);
          materialize(chosen);
        }
      }
      return;
    }

    const ui pSize = popcount(p);
    ui maximumSize = 0;
    if (!cliqueSizeWithExtension(cliqueSize, pSize, maximumSize))
      return;
    if (maximumSize < minCliqueSize)
      return;

    ui pivot = 0;
    ui bestScore = 0;
    bool havePivot = false;
    bool pIsClique = true;
    bool complementIsMatching = true;
    ui complementDegreeSum = 0;

    ull candidates = p | x;
    while (candidates != 0) {
      const ui index = static_cast<ui>(__builtin_ctzll(candidates));
      const ull bit = 1ULL << index;
      candidates &= candidates - 1;

      const ui score = popcount(p & neighbors[index]);
      if (!havePivot || score > bestScore) {
        pivot = index;
        bestScore = score;
        havePivot = true;
      }

      // An original or recursively processed X vertex adjacent to all of P
      // blocks every continuation in this BK state.
      if ((x & bit) != 0 && score == pSize)
        return;

      if ((p & bit) != 0) {
        const ui complementDegree = pSize - 1 - score;
        complementDegreeSum += complementDegree;
        pIsClique = pIsClique && complementDegree == 0;
        complementIsMatching =
            complementIsMatching && complementDegree <= 1;
      }
    }

    // P is complete and no X vertex covers it, so R union P is its sole
    // maximal continuation.
    if (pIsClique) {
      if (addCliques(1, maximumSize)) {
        saveWitness(chosen | p);
        materialize(chosen | p);
      }
      return;
    }

    if (x == 0 && complementIsMatching) {
      const ui missingEdges = complementDegreeSum >> 1;
      const ui extensionSize = pSize - missingEdges;
      ui maximalSize = 0;
      if (!cliqueSizeWithExtension(cliqueSize, extensionSize, maximalSize))
        return;
      if (maximalSize < minCliqueSize)
        return;
      if (missingEdges >= 64) {
        overflow = true;
        return;
      }

      const ull count = 1ULL << missingEdges;
      if (addCliques(count, maximalSize)) {
        saveWitness(chosen | matchingWitness(p));
        if (cliqueSink != nullptr)
          materializeMatching(chosen, p);
      }
      return;
    }

    ull branch = p & ~neighbors[pivot];
    while (branch != 0 && !overflow) {
      const ui index = static_cast<ui>(__builtin_ctzll(branch));
      const ull bit = 1ULL << index;
      branch &= branch - 1;

      enumerate(p & neighbors[index], x & neighbors[index], cliqueSize + 1,
                chosen | bit);
      p &= ~bit;
      x |= bit;
    }
  }

public:
  OneWordSubtreeSolver(const FastAdjacencyHash &adjacency,
                       const std::vector<ui> &localVertices,
                       const std::vector<ui> *prefix, ui outputThreshold,
                       const FastCliqueSink *outputSink)
      : vertices(localVertices), cliquePrefix(prefix), cliqueSink(outputSink),
        minCliqueSize(std::max<ui>(1, outputThreshold)) {
    for (ui i = 0; i < vertices.size(); ++i) {
      ull mask = 0;
      for (ui j = 0; j < vertices.size(); ++j) {
        if (i != j && adjacency.contains(vertices[i], vertices[j]))
          mask |= 1ULL << j;
      }
      neighbors[i] = mask;
    }
  }

  FastLocalBitsetResult run(ui pSize, ui cliqueSize) {
    const ull p = pSize == 64 ? std::numeric_limits<ull>::max()
                              : ((1ULL << pSize) - 1);
    const ui xSize = static_cast<ui>(vertices.size()) - pSize;
    const ull x = xSize == 0
                      ? 0
                      : (xSize == 64
                             ? std::numeric_limits<ull>::max()
                             : (((1ULL << xSize) - 1) << pSize));

    enumerate(p, x, cliqueSize, 0);
    result.handled = !overflow;
    if (overflow) {
      // No partial aggregate may escape: handled=false means the caller will
      // enumerate this subtree from its original, read-only P/X state.
      result.found = false;
      result.cliqueCount = 0;
      result.maxCliqueSize = 0;
      result.witness.clear();
      pendingCliques.clear();
    } else if (cliqueSink != nullptr) {
      for (const std::vector<ui> &clique : pendingCliques)
        (*cliqueSink)(clique);
    }
    return result;
  }
};

} // namespace

FastLocalBitsetResult solveFastLocalBitsetSubtree(
    const FastAdjacencyHash &adjacency, const std::vector<ui> &p,
    const std::vector<ui> &x, ui cliqueSize,
    const std::vector<ui> *cliquePrefix, ui minCliqueSize,
    const FastCliqueSink *cliqueSink) {
  FastLocalBitsetResult declined;
  if (p.size() > 64 || x.size() > 64 - p.size())
    return declined;

  std::vector<ui> vertices;
  vertices.reserve(p.size() + x.size());
  vertices.insert(vertices.end(), p.begin(), p.end());
  vertices.insert(vertices.end(), x.begin(), x.end());

  OneWordSubtreeSolver solver(adjacency, vertices, cliquePrefix,
                              minCliqueSize, cliqueSink);
  return solver.run(static_cast<ui>(p.size()), cliqueSize);
}
