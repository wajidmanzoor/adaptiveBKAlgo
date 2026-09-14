#pragma once

#include "checked_count.h"
#include "common.h"
#include "graph.h"

#include <memory>

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

  struct CcrBitState {
    vector<ull> core;
    vector<ull> residual;
    vector<ull> excludedCandidates;
    vector<ui> excludedExternal;
    vector<ui> branchRoots;
    ui coreCount = 0;
    ui residualCount = 0;
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
  ull solverWorkBudget;
  bool solverWorkBudgetEnabled;
  ui minCliqueSize;
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
  // Reusable local CCRMCE representation. Every row in ccrNeighborBits is a
  // bitset over the current Q. Rows [0,Q.size()) belong to Q and the
  // remaining rows belong to the initial external X when materialized.
  // Keeping X sparse while C/P/X-within-Q are bitsets avoids repeated
  // allocation/copy costs without materializing a full (Q union X)^2 graph.
  vector<CcrBitState> ccrStates;
  vector<ui> ccrQVertices;
  vector<ui> ccrXVertices;
  vector<ull> ccrNeighborBits;
  size_t ccrWordCount = 0;
  bool ccrExternalRowsMaterialized = false;
  vector<ui> ccrHeap;
  vector<ui> ccrHeapPosition;
  vector<ui> ccrDegree;
  // Reusable vertex-to-local-Q map for local induced-graph construction. It
  // is separate from the seed solver's E-index map; E's map is reused for X.
  vector<ui> ccrIndex;
  vector<ui> ccrIndexStamp;
  ui ccrIndexToken;
  // Reused across top-level CCRMCE invocations. Calls are sequential (the
  // recursive search only mutates ccrStates), so retaining these capacities
  // removes repeated branch-local allocation without changing search state.
  vector<ui> ccrCommon;
  vector<ui> ccrCommonScratch;
  vector<ui> ccrBranchX;
  vector<ui> ccrBranchR;
  vector<ui> coveringCliqueScratch;
  ull ccrFindOneCalls;
  ull ccrFullCalls;
  ull ccrFindOneStates;
  ull ccrFullStates;
  ull ccrCoreExtractions;
  ull ccrCoreVertices;
  ull ccrResidualVertices;
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

  void commonExpandInto(vector<ui> &out, const vector<ui> &E,
                        const vector<ui> &S);
  vector<vector<ui>> buildHitSets(const vector<ui> &E,
                                  const vector<ui> &cliqueIds,
                                  ui maxHitSets = UINT_MAX);
  vector<vector<ui>> singletonBranches(const vector<ui> &E);

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

  bool collectAllCoveringCliques(const vector<ui> &M,
                                 vector<ui> &result);
  vector<vector<ui>>
  generateExactSiblingSets(const vector<ui> &E,
                           const vector<ui> &coveringCliqueIds,
                           bool *usePivotFallback = nullptr);
  bool findOnePure(const vector<ui> &M, const vector<ui> &Q,
                   vector<ui> &found);
  bool findOnePureRecursive(vector<ui> &R, CcrBitState &state,
                            bool fromP, vector<ui> &found, size_t depth);
  bool runSmallCcrBranch(const vector<ui> &M, const vector<ui> &Q,
                         const vector<ui> &common, vector<ui> *found);
  void prepareCcrBranch(const vector<ui> &Q, const vector<ui> &X,
                        bool materializeExternalRows);
  void normalizeCcrState(vector<ui> &R, CcrBitState &state,
                         bool &fromP);
  void pruneCcrExcludedWithoutResidualNeighbors(CcrBitState &state);
  void buildCcrChildState(CcrBitState &child,
                          const CcrBitState &parent, ui local) const;
  bool ccrCoreUnionIsMaximal(const CcrBitState &state) const;
  size_t selectCcrPivot(const CcrBitState &state) const;
  void buildCcrBranchRoots(CcrBitState &state, size_t pivot) const;
  void enumerateAllPureBranch(const vector<ui> &M, const vector<ui> &Q);
  void enumerateAllPureBranchRecursive(vector<ui> &R, CcrBitState &state,
                                       bool fromP, size_t depth);
  bool recordPureClique(vector<ui> C);

public:
  explicit ReorderSib(Graph &g, ui minCliqueSize = 3);
  void findAllMaximalCliquesPure();
  void setSolverWorkBudget(ull budget) {
    solverWorkBudget = budget;
    solverWorkBudgetEnabled = true;
  }
  void clearSolverWorkBudget() {
    solverWorkBudget = 0;
    solverWorkBudgetEnabled = false;
  }
  vector<vector<ui>> getCliques() const;
};
