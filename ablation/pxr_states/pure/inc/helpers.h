#pragma once

#include "checked_count.h"
#include "common.h"
#include "graph.h"

#include <memory>

enum class SibMethod {
  BACKTRACKING,
  OPTIMIZED
};

class ReorderSib {
private:
  struct AdjacencyRow {
    const ui *data;
    size_t length;

    const ui *begin() const { return data; }
    const ui *end() const { return data + length; }
    size_t size() const { return length; }
    bool empty() const { return length == 0; }
    ui operator[](size_t index) const { return data[index]; }
  };

  struct PureBranch {
    vector<ui> mustin;
    vector<ui> expandTo;
  };

  ui n;
  // Reordered graph in CSR form: row u is
  // adjVertices[adjOffsets[u]..adjOffsets[u + 1]).
  vector<ui> adjVertices;
  vector<size_t> adjOffsets;
  vector<ui> firstForwardNeighbor;
  // Low-degree rows use binary search in CSR. Only high-degree rows pay
  // for a hash table, avoiding one heavyweight unordered_set per vertex.
  vector<unique_ptr<unordered_set<ui>>> adjHash;
  ull cliqueCount;
  ull dupBlocked;
  size_t maxCliqueSize;
  ull checksCount;
  ull findOnePxrStates;
  ull fullPxrStates;
  ull et1EnumeratedStates;
  ull et2EnumeratedStates;
  ull et3EnumeratedStates;
  ull solverWorkBudget;
  bool solverWorkBudgetEnabled;
  ull solverBudgetFallbacks;
  ull solverCertifiedBudgetFallbacks;
  ull findOne2PlexTerminals;
  ull findOne3PlexTerminals;
  ull fullPxr2PlexTerminals;
  ull fullPxr3PlexTerminals;
  ull worklistPushes;
  ull worklistPops;
  size_t maximumWorklistSize;
  ull minSizePrunedBranches;
  ull coverLookupCalls;
  ull findOneCalls;
  ull findOneSuccesses;
  ull findOneCliqueSizeTotal;
  size_t maximumFindOneCliqueSize;
  ull seedSolverCalls;
  ull generatedBranches;
  ull fullPxrFallbackBranches;
  ull smallQFullPxrBranches;
  array<ull, 8> poppedQSizeBuckets;
  SibMethod method;
  ui minCliqueSize;
  ui smallQFullPxrThreshold;
  // Compact append-only clique arena. Clique i occupies
  // cliqueVertices[cliqueOffsets[i]..cliqueOffsets[i + 1]).
  vector<ui> cliqueVertices;
  vector<size_t> cliqueOffsets{0};
  // The search and cover index use internal permuted labels. External
  // validation restores original graph IDs through this inverse permutation.
  vector<ui> internalToOriginal;
  // Pure stores every emitted clique in one per-vertex posting list. The old
  // representation retained a redundant level dimension whose only live
  // bucket was always level zero.
  vector<vector<ui>> cliqueIdsByVertex;
  // Open-addressed hash slots point to exact same-hash chains. Distinct
  // cliques with a colliding 64-bit hash remain separate and are resolved by
  // comparing their complete arena contents.
  vector<ull> emittedHashKeys;
  vector<size_t> emittedHashHeads;
  vector<size_t> emittedHashNext;
  size_t emittedHashSlotsUsed = 0;

  vector<ui> eIndex;             // reusable vertex -> local E index map
  vector<ui> eIndexStamp;        // stamp for entries valid in current solver call
  ui eIndexToken;
  // P/X child buffers are indexed by recursion depth and reused across
  // sibling calls and top-level branches.
  vector<vector<ui>> pxrPBuffers;
  vector<vector<ui>> pxrXBuffers;
  // Reused by commonExpandInto so sibling generation does not allocate an
  // ordering and intersection scratch vector for every branch.
  vector<ui> commonExpandOrder;
  vector<ui> commonExpandScratch;

  AdjacencyRow adjacentVertices(ui u) const {
    static constexpr ui emptyRow = 0;
    const ui *base = adjVertices.empty() ? &emptyRow : adjVertices.data();
    return {base + adjOffsets[u], adjOffsets[u + 1] - adjOffsets[u]};
  }

  void intersectInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void intersectInto(vector<ui> &out, const vector<ui> &A, AdjacencyRow B);
  void intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                              const vector<ui> &B,
                              const vector<ui> &exclude);
  void intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                              AdjacencyRow B, const vector<ui> &exclude);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  vector<ui> setDiffStoredClique(const vector<ui> &A, ui cliqueId) const;
  void setDiffInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void setDiffInto(vector<ui> &out, const vector<ui> &A, AdjacencyRow B);
  void unionSetInto(vector<ui> &out, const vector<ui> &A,
                    const vector<ui> &B);
  size_t storedCliqueCount() const { return cliqueOffsets.size() - 1; }
  size_t storedCliqueSize(ui cliqueId) const {
    return cliqueOffsets[cliqueId + 1] - cliqueOffsets[cliqueId];
  }
  bool storedCliqueContains(ui cliqueId, const vector<ui> &subset) const;
  bool storedCliqueEquals(ui cliqueId, const vector<ui> &clique) const;
  void rehashEmittedCliqueIndex(size_t capacity);
  size_t emittedCliqueHashSlot(ull hash) const;

  bool hitsAll(const vector<ui> &S, const vector<vector<ui>> &hitSets);
  void commonExpandInto(vector<ui> &out, const vector<ui> &E,
                        const vector<ui> &S);
  vector<vector<ui>> buildHitSets(const vector<ui> &E,
                                  const vector<ui> &cliqueIds,
                                  ui maxHitSets = UINT_MAX);
  vector<vector<ui>> singletonBranches(const vector<ui> &E);

  vector<vector<ui>> backtrackingBranchBound(const vector<ui> &E,
                                             const vector<vector<ui>> &hitSets);
  vector<vector<ui>> efficientHittingSet(const vector<ui> &E,
                                         vector<vector<ui>> hitSets,
                                         bool *usePivotFallback = nullptr);
  vector<vector<ui>> efficientHittingSetDirect(
      const vector<ui> &E, const vector<ui> &coveringCliqueIds,
      bool *usePivotFallback = nullptr);
  bool adj(ui u, ui v) const {
    const auto &hash = adjHash[u];
    if (hash)
      return hash->find(v) != hash->end();
    const AdjacencyRow row = adjacentVertices(u);
    return binary_search(row.begin(), row.end(), v);
  }

  vector<ui> collectAllCoveringCliques(const vector<ui> &M);
  vector<vector<ui>>
  generateExactSiblingSets(const vector<ui> &E,
                           const vector<ui> &coveringCliqueIds,
                           bool *usePivotFallback = nullptr);
  bool findOnePure(const vector<ui> &M, const vector<ui> &Q,
                   vector<ui> &found);
  bool findOnePureRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X,
                            vector<ui> &found, size_t depth);
  ui pureNeighborsInP(ui u, const vector<ui> &P) const;
  void scanPurePXRState(const vector<ui> &P, const vector<ui> &X,
                        ui &pivot, ui &minPScore, ui &universalP,
                        bool &xUniversal);
  void pureMatchingParts(const vector<ui> &P, vector<ui> &forced,
                         vector<pair<ui, ui>> &missingEdges) const;
  void enumerateSmallPureBranch(const vector<ui> &M,
                                const vector<ui> &Q);
  void enumerateAllPureBranch(const vector<ui> &M, const vector<ui> &Q);
  void enumerateAllPureBranchRecursive(vector<ui> &R, vector<ui> &P,
                                       vector<ui> &X, size_t depth);
  bool recordPureClique(vector<ui> C);

public:
  ReorderSib(Graph &g, SibMethod method = SibMethod::OPTIMIZED,
             ui minCliqueSize = 3);
  void findAllMaximalCliquesPure();
  void setSolverWorkBudget(ull budget) {
    solverWorkBudget = budget;
    solverWorkBudgetEnabled = true;
  }
  void clearSolverWorkBudget() {
    solverWorkBudget = 0;
    solverWorkBudgetEnabled = false;
  }
  void setSmallQFullPxrThreshold(ui threshold) {
    smallQFullPxrThreshold = threshold;
  }
  vector<vector<ui>> getCliques() const;
};
