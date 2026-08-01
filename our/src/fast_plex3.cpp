#include "../inc/fast_plex3.h"
#include "../inc/checked_count.h"

#include <array>
#include <cstddef>
#include <limits>

namespace {

constexpr size_t SMALL_SIZE_LIMIT = 2;

struct Counts {
  ull total = 0;
  std::array<ull, SMALL_SIZE_LIMIT + 1> bySize{};
};

struct ComponentStats {
  Counts counts;
  ui maximumSize = 0;
};

Counts singletonSequence(ui selected) {
  Counts result;
  result.total = 1;
  if (selected <= SMALL_SIZE_LIMIT)
    result.bySize[selected] = 1;
  return result;
}

bool appendCounts(Counts &destination, const Counts &source,
                  ui selectedIncrement) {
  ull sum = 0;
  if (!tryAddUll(destination.total, source.total, sum))
    return false;
  destination.total = sum;

  for (size_t selected = 0;
       selected + selectedIncrement <= SMALL_SIZE_LIMIT; ++selected) {
    const size_t target = selected + selectedIncrement;
    if (!tryAddUll(destination.bySize[target], source.bySize[selected], sum))
      return false;
    destination.bySize[target] = sum;
  }
  return true;
}

// Count binary strings that encode maximal independent sets.  A selected
// vertex is 1.  Independence forbids 11; maximality forbids an undominated
// zero, i.e. 000 internally and 00 at either endpoint of a path.
bool countPath(size_t length, ComponentStats &result) {
  if (length == 0)
    return false;
  if (length == 1) {
    result.counts = singletonSequence(1);
    result.maximumSize = 1;
    return true;
  }

  // The low bit is the newest bit.  Only 01 and 10 satisfy both endpoint
  // domination and independence for the first two path vertices.
  std::array<Counts, 4> states{};
  states[1] = singletonSequence(1);
  states[2] = singletonSequence(1);

  for (size_t at = 2; at < length; ++at) {
    std::array<Counts, 4> next{};
    for (ui state = 0; state < states.size(); ++state) {
      if (states[state].total == 0)
        continue;
      const ui previousPrevious = state >> 1;
      const ui previous = state & 1U;
      for (ui selected = 0; selected <= 1; ++selected) {
        if (previous != 0 && selected != 0)
          continue;
        if (previous == 0 && previousPrevious == 0 && selected == 0)
          continue;
        const ui nextState = (previous << 1) | selected;
        if (!appendCounts(next[nextState], states[state], selected))
          return false;
      }
    }
    states = next;
  }

  // The final endpoint must be selected or dominated by its predecessor.
  if (!appendCounts(result.counts, states[1], 0) ||
      !appendCounts(result.counts, states[2], 0))
    return false;
  result.maximumSize = static_cast<ui>((length + 1) / 2);
  return true;
}

bool countCycle(size_t length, ComponentStats &result) {
  if (length < 3)
    return false;

  // State = first two bits followed by the current last two bits.  Retaining
  // the first pair lets the final transition validate the wraparound edge
  // and domination of both boundary vertices.
  std::array<Counts, 16> states{};
  for (ui firstPair = 0; firstPair < 3; ++firstPair) {
    const ui selected = ((firstPair >> 1) & 1U) + (firstPair & 1U);
    states[firstPair * 4 + firstPair] = singletonSequence(selected);
  }

  for (size_t at = 2; at < length; ++at) {
    std::array<Counts, 16> next{};
    for (ui state = 0; state < states.size(); ++state) {
      if (states[state].total == 0)
        continue;
      const ui firstPair = state >> 2;
      const ui lastPair = state & 3U;
      const ui previousPrevious = lastPair >> 1;
      const ui previous = lastPair & 1U;
      for (ui selected = 0; selected <= 1; ++selected) {
        if (previous != 0 && selected != 0)
          continue;
        if (previous == 0 && previousPrevious == 0 && selected == 0)
          continue;
        const ui nextPair = (previous << 1) | selected;
        if (!appendCounts(next[firstPair * 4 + nextPair], states[state],
                          selected))
          return false;
      }
    }
    states = next;
  }

  for (ui state = 0; state < states.size(); ++state) {
    const Counts &source = states[state];
    if (source.total == 0)
      continue;

    const ui firstPair = state >> 2;
    const ui lastPair = state & 3U;
    const ui first = firstPair >> 1;
    const ui second = firstPair & 1U;
    const ui previous = lastPair >> 1;
    const ui last = lastPair & 1U;

    if (first != 0 && last != 0)
      continue;
    if (last == 0 && previous == 0 && first == 0)
      continue;
    if (first == 0 && last == 0 && second == 0)
      continue;
    if (!appendCounts(result.counts, source, 0))
      return false;
  }

  result.maximumSize = static_cast<ui>(length / 2);
  return true;
}

bool combineSmallCounts(
    const std::array<ull, SMALL_SIZE_LIMIT + 1> &left,
    const std::array<ull, SMALL_SIZE_LIMIT + 1> &right,
    std::array<ull, SMALL_SIZE_LIMIT + 1> &result) {
  result.fill(0);
  for (size_t i = 0; i <= SMALL_SIZE_LIMIT; ++i) {
    for (size_t j = 0; i + j <= SMALL_SIZE_LIMIT; ++j) {
      ull product = 0;
      ull sum = 0;
      if (!tryMultiplyUll(left[i], right[j], product) ||
          !tryAddUll(result[i + j], product, sum))
        return false;
      result[i + j] = sum;
    }
  }
  return true;
}

void enumerateComponentSelections(
    const std::vector<size_t> &order, bool cycle,
    const std::vector<ui> &vertices,
    std::vector<std::vector<ui>> &selections) {
  std::vector<unsigned char> selected(order.size(), 0);
  std::function<void(size_t)> enumerate = [&](size_t at) {
    if (at != order.size()) {
      selected[at] = 0;
      enumerate(at + 1);
      if (at == 0 || selected[at - 1] == 0) {
        selected[at] = 1;
        enumerate(at + 1);
        selected[at] = 0;
      }
      return;
    }

    if (cycle && selected.front() != 0 && selected.back() != 0)
      return;
    for (size_t i = 0; i < selected.size(); ++i) {
      if (selected[i] != 0)
        continue;
      const bool left =
          i != 0 ? selected[i - 1] != 0
                 : cycle && selected.back() != 0;
      const bool right =
          i + 1 != selected.size()
              ? selected[i + 1] != 0
              : cycle && selected.front() != 0;
      if (!left && !right)
        return;
    }

    std::vector<ui> extension;
    for (size_t i = 0; i < selected.size(); ++i)
      if (selected[i] != 0)
        extension.push_back(vertices[order[i]]);
    selections.push_back(std::move(extension));
  };
  enumerate(0);
}

template <typename Contains>
FastPlex3Result solveFastPlex3SubtreeImpl(
    Contains contains, const std::vector<ui> &p,
    ui cliqueSize, const std::vector<ui> *cliquePrefix,
    ui minCliqueSize, const FastCliqueSink *cliqueSink) {
  FastPlex3Result result;
  minCliqueSize = std::max<ui>(1, minCliqueSize);
  // The dynamic program retains exact counts for candidate extensions of
  // sizes 0, 1, and 2 so the default tau=3 filter is constant-space.  Larger
  // output thresholds fall back to ordinary BK rather than changing either
  // the complement-degree-two recognition rule or its hot-path footprint.
  if (minCliqueSize > SMALL_SIZE_LIMIT + 1)
    return result;
  const size_t pSize = p.size();
  if (pSize > std::numeric_limits<ui>::max())
    return result;

  // Store the at-most-two complement neighbors of each P vertex.  Finding a
  // third proves this is not a 3-plex, so the caller keeps ordinary BK.
  const size_t noVertex = pSize;
  std::vector<std::array<size_t, 2>> complement(pSize);
  std::vector<unsigned char> complementDegree(pSize, 0);
  for (size_t i = 0; i < pSize; ++i) {
    for (size_t j = i + 1; j < pSize; ++j) {
      ++result.checksCount;
      if (contains(p[i], p[j]))
        continue;
      if (complementDegree[i] == 2 || complementDegree[j] == 2)
        return result;
      complement[i][complementDegree[i]++] = j;
      complement[j][complementDegree[j]++] = i;
    }
  }

  ull totalCount = 1;
  std::array<ull, SMALL_SIZE_LIMIT + 1> smallCounts{{1, 0, 0}};
  ui maximumCandidateSize = 0;
  std::vector<ui> maximumCandidates;
  maximumCandidates.reserve(pSize);

  std::vector<unsigned char> visited(pSize, 0);
  std::vector<size_t> stack;
  std::vector<size_t> component;
  std::vector<size_t> order;
  std::vector<std::vector<std::vector<ui>>> componentSelections;
  stack.reserve(pSize);
  component.reserve(pSize);
  order.reserve(pSize);

  for (size_t seed = 0; seed < pSize; ++seed) {
    if (visited[seed] != 0)
      continue;

    stack.clear();
    component.clear();
    stack.push_back(seed);
    visited[seed] = 1;
    while (!stack.empty()) {
      const size_t vertex = stack.back();
      stack.pop_back();
      component.push_back(vertex);
      for (ui at = 0; at < complementDegree[vertex]; ++at) {
        const size_t neighbor = complement[vertex][at];
        if (visited[neighbor] == 0) {
          visited[neighbor] = 1;
          stack.push_back(neighbor);
        }
      }
    }

    bool cycle = true;
    size_t start = component.front();
    for (size_t vertex : component) {
      if (complementDegree[vertex] < 2) {
        cycle = false;
        if (complementDegree[vertex] <= 1)
          start = vertex;
      }
    }

    order.clear();
    size_t previous = noVertex;
    size_t current = start;
    for (;;) {
      order.push_back(current);
      size_t next = noVertex;
      for (ui at = 0; at < complementDegree[current]; ++at) {
        const size_t candidate = complement[current][at];
        if (candidate != previous) {
          next = candidate;
          break;
        }
      }
      if (next == noVertex || (cycle && next == start))
        break;
      previous = current;
      current = next;
    }
    if (order.size() != component.size())
      return result;

    ComponentStats stats;
    if (!(cycle ? countCycle(order.size(), stats)
                : countPath(order.size(), stats)))
      return result;
    if (cliqueSink != nullptr) {
      componentSelections.emplace_back();
      enumerateComponentSelections(order, cycle, p,
                                   componentSelections.back());
      if (componentSelections.back().empty())
        return result;
    }

    ull combinedTotal = 0;
    if (!tryMultiplyUll(totalCount, stats.counts.total, combinedTotal))
      return result;
    totalCount = combinedTotal;

    std::array<ull, SMALL_SIZE_LIMIT + 1> combinedSmall{};
    if (!combineSmallCounts(smallCounts, stats.counts.bySize,
                            combinedSmall))
      return result;
    smallCounts = combinedSmall;

    if (stats.maximumSize >
        std::numeric_limits<ui>::max() - maximumCandidateSize)
      return result;
    maximumCandidateSize += stats.maximumSize;

    size_t witnessLimit = order.size();
    if (cycle && (order.size() & 1U) != 0)
      --witnessLimit;
    for (size_t at = 0; at < witnessLimit; at += 2)
      maximumCandidates.push_back(p[order[at]]);
  }

  if (maximumCandidateSize >
      std::numeric_limits<ui>::max() - cliqueSize)
    return result;
  const ui maximumCliqueSize = cliqueSize + maximumCandidateSize;

  ull excluded = 0;
  if (cliqueSize < minCliqueSize) {
    const ui maximumExcludedCandidateSize =
        minCliqueSize - 1 - cliqueSize;
    for (ui selected = 0; selected <= maximumExcludedCandidateSize;
         ++selected) {
      ull sum = 0;
      if (!tryAddUll(excluded, smallCounts[selected], sum))
        return result;
      excluded = sum;
    }
  }
  if (excluded > totalCount)
    return result;

  result.handled = true;
  result.cliqueCount = totalCount - excluded;
  result.found = result.cliqueCount != 0;
  if (!result.found)
    return result;

  result.maxCliqueSize = maximumCliqueSize;
  if (cliquePrefix != nullptr)
    result.witness = *cliquePrefix;
  result.witness.insert(result.witness.end(), maximumCandidates.begin(),
                        maximumCandidates.end());

  if (cliqueSink != nullptr) {
    std::vector<std::vector<ui>> pendingCliques;
    std::vector<ui> extension;
    std::function<void(size_t)> combine = [&](size_t componentIndex) {
      if (componentIndex == componentSelections.size()) {
        if (cliqueSize + extension.size() < minCliqueSize)
          return;
        std::vector<ui> clique;
        if (cliquePrefix != nullptr)
          clique = *cliquePrefix;
        clique.insert(clique.end(), extension.begin(), extension.end());
        std::sort(clique.begin(), clique.end());
        pendingCliques.push_back(std::move(clique));
        return;
      }

      for (const std::vector<ui> &selection :
           componentSelections[componentIndex]) {
        const size_t oldSize = extension.size();
        extension.insert(extension.end(), selection.begin(), selection.end());
        combine(componentIndex + 1);
        extension.resize(oldSize);
      }
    };
    combine(0);
    if (pendingCliques.size() != result.cliqueCount)
      return FastPlex3Result{};
    for (const std::vector<ui> &clique : pendingCliques)
      (*cliqueSink)(clique);
  }
  return result;
}

} // namespace

FastPlex3Result solveFastPlex3Subtree(
    const FastAdjacencyHash &adjacency, const std::vector<ui> &p,
    ui cliqueSize, const std::vector<ui> *cliquePrefix,
    ui minCliqueSize, const FastCliqueSink *cliqueSink) {
  return solveFastPlex3SubtreeImpl(
      [&](ui u, ui v) { return adjacency.contains(u, v); }, p, cliqueSize,
      cliquePrefix, minCliqueSize, cliqueSink);
}

FastPlex3Result solveFastPlex3Subtree(
    const std::vector<std::unordered_set<ui>> &adjacency,
    const std::vector<ui> &p, ui cliqueSize,
    const std::vector<ui> *cliquePrefix, ui minCliqueSize,
    const FastCliqueSink *cliqueSink) {
  return solveFastPlex3SubtreeImpl(
      [&](ui u, ui v) { return adjacency[u].count(v) != 0; }, p, cliqueSize,
      cliquePrefix, minCliqueSize, cliqueSink);
}
