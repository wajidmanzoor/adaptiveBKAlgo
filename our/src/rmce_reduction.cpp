#include "../inc/rmce_reduction.h"

#include <limits>
#include <queue>
#include <set>
#include <stdexcept>

namespace {

void incrementOrThrow(ull &value, const char *what) {
  if (value == numeric_limits<ull>::max())
    throw overflow_error(what);
  ++value;
}

} // namespace

RmceReductionResult applyRmceReduction(const Graph &graph, ui minCliqueSize,
                                       bool collectCliqueIdentities) {
  const ui n = graph.n;
  minCliqueSize = max<ui>(1, minCliqueSize);

  vector<set<ui>> adjacency(n);
  vector<ui> originalDegree(n, 0);
  vector<unsigned char> active(n, 1);
  for (ui v = 0; v < n; ++v) {
    adjacency[v].insert(graph.neighbors.begin() + graph.offset[v],
                        graph.neighbors.begin() + graph.offset[v + 1]);
    adjacency[v].erase(v);
    originalDegree[v] = static_cast<ui>(adjacency[v].size());
  }

  RmceReductionResult result;
  set<vector<ui>> emitted;
  auto emit = [&](vector<ui> clique) {
    sort(clique.begin(), clique.end());
    if (clique.size() < minCliqueSize) {
      incrementOrThrow(result.counters.discardedUnderThreshold,
                       "RMCE under-threshold counter exceeds uint64_t");
      return;
    }
    if (collectCliqueIdentities) {
      if (!emitted.insert(clique).second) {
        incrementOrThrow(result.counters.duplicateDirectOutputs,
                         "RMCE duplicate counter exceeds uint64_t");
        return;
      }
      result.directlyEmittedCliques.push_back(clique);
    }
    incrementOrThrow(result.directlyEmittedCount,
                     "RMCE direct output count exceeds uint64_t");
    incrementOrThrow(result.counters.directlyEmitted,
                     "RMCE direct output counter exceeds uint64_t");
    result.maximumCliqueSize = max(result.maximumCliqueSize, clique.size());
  };

  queue<ui> lowDegree;
  vector<unsigned char> queued(n, 0);
  auto enqueueIfLow = [&](ui v) {
    if (v < n && active[v] != 0 && adjacency[v].size() <= 2 &&
        queued[v] == 0) {
      queued[v] = 1;
      lowDegree.push(v);
    }
  };
  for (ui v = 0; v < n; ++v)
    enqueueIfLow(v);

  auto removeEdge = [&](ui u, ui v) {
    adjacency[u].erase(v);
    adjacency[v].erase(u);
    enqueueIfLow(u);
    enqueueIfLow(v);
  };

  auto commonNeighborCount = [&](ui u, ui v) {
    const set<ui> &a = adjacency[u];
    const set<ui> &b = adjacency[v];
    auto i = a.begin();
    auto j = b.begin();
    ui count = 0;
    while (i != a.end() && j != b.end()) {
      if (*i < *j) {
        ++i;
      } else if (*j < *i) {
        ++j;
      } else {
        if (active[*i] != 0)
          ++count;
        ++i;
        ++j;
      }
    }
    return count;
  };

  auto processLowDegree = [&]() {
    bool changed = false;
    while (!lowDegree.empty()) {
      const ui u = lowDegree.front();
      lowDegree.pop();
      queued[u] = 0;
      if (active[u] == 0 || adjacency[u].size() > 2)
        continue;

      changed = true;
      const vector<ui> neighbors(adjacency[u].begin(), adjacency[u].end());
      if (neighbors.empty()) {
        incrementOrThrow(result.counters.degree0Vertices,
                         "RMCE degree-zero counter exceeds uint64_t");
        // Vertices isolated by earlier reductions are silent because their
        // incident maximal cliques have already been emitted.
        if (originalDegree[u] == 0)
          emit({u});
        active[u] = 0;
        continue;
      }

      if (neighbors.size() == 1) {
        incrementOrThrow(result.counters.degree1Vertices,
                         "RMCE degree-one counter exceeds uint64_t");
        const ui v = neighbors[0];
        emit({u, v});
        removeEdge(u, v);
        active[u] = 0;
        adjacency[u].clear();
        enqueueIfLow(v);
        continue;
      }

      incrementOrThrow(result.counters.degree2Vertices,
                       "RMCE degree-two counter exceeds uint64_t");
      const ui v = neighbors[0];
      const ui w = neighbors[1];
      const bool vwEdge = adjacency[v].find(w) != adjacency[v].end();
      if (!vwEdge) {
        emit({u, v});
        emit({u, w});
      } else {
        emit({u, v, w});
        if (commonNeighborCount(v, w) == 1)
          removeEdge(v, w);
      }
      removeEdge(u, v);
      removeEdge(u, w);
      active[u] = 0;
      adjacency[u].clear();
      enqueueIfLow(v);
      enqueueIfLow(w);
    }
    return changed;
  };

  // Removing a non-triangle edge can expose another low-degree vertex, so
  // alternate both global RMCE phases until neither changes the graph.
  while (true) {
    bool changed = processLowDegree();
    vector<pair<ui, ui>> edges;
    for (ui u = 0; u < n; ++u) {
      if (active[u] == 0)
        continue;
      for (ui v : adjacency[u])
        if (u < v && active[v] != 0)
          edges.emplace_back(u, v);
    }
    for (const auto &[u, v] : edges) {
      if (active[u] == 0 || active[v] == 0 ||
          adjacency[u].find(v) == adjacency[u].end())
        continue;
      if (commonNeighborCount(u, v) == 0) {
        emit({u, v});
        incrementOrThrow(result.counters.nontriangleEdges,
                         "RMCE non-triangle counter exceeds uint64_t");
        removeEdge(u, v);
        changed = true;
      }
    }
    if (!changed && lowDegree.empty())
      break;
  }

  vector<ui> oldToNew(n, UINT_MAX);
  for (ui old = 0; old < n; ++old) {
    if (active[old] != 0 && !adjacency[old].empty()) {
      oldToNew[old] = static_cast<ui>(result.residualToOriginal.size());
      result.residualToOriginal.push_back(old);
    }
  }

  vector<pair<ui, ui>> residualEdges;
  for (ui oldU = 0; oldU < n; ++oldU) {
    const ui u = oldToNew[oldU];
    if (u == UINT_MAX)
      continue;
    for (ui oldV : adjacency[oldU]) {
      const ui v = oldToNew[oldV];
      if (v != UINT_MAX && u < v)
        residualEdges.emplace_back(u, v);
    }
  }
  result.graph =
      Graph(static_cast<ui>(result.residualToOriginal.size()), residualEdges);
  return result;
}
