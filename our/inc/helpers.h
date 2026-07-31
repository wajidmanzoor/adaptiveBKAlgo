#pragma once

#include "checked_count.h"
#include "common.h"
#include "graph.h"

enum class DegOrder { ORIGINAL, ASCENDING, DESCENDING };

struct ReorderSibTestAccess;

// Baseline adjacency-list Bron--Kerbosch with pivoting.
class PivotBK {
private:
  ui n;
  vector<vector<ui>> adjList;
  ull cliqueCount;
  ui maxCliqueSize;
  ull checksCount;

  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);

public:
  explicit PivotBK(Graph &g, DegOrder order = DegOrder::ASCENDING);

  void findAllMaximalCliques();
  ull getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
  ull getChecksCount() const { return checksCount; }
};

// Exact Pure-ReorderSib implementation.
// 128-constraint masks, PXR/ET terminal handling, and fail-closed budget
// fallback to exhaustive pivot enumeration.
class ReorderSib {
private:
  friend struct ReorderSibTestAccess;

  struct PureBranch {
    vector<ui> mustin;
    vector<ui> expandTo;
  };

  ui n;
  vector<vector<ui>> adjList;
  vector<vector<ui>> adjList2;
  vector<unordered_set<ui>> adjSet;
  ull cliqueCount;
  ull dupBlocked;
  size_t maxCliqueSize;
  ull checksCount;
  ull solverWorkBudget;
  ull solverBudgetFallbacks;
  ui minCliqueSize;

  vector<vector<ui>> allCliques;
  vector<ui> internalToOriginal;
  vector<vector<vector<ui>>> cliquesByVertexByLevel;
  vector<ull> cliqueCountByVertex;
  unordered_set<string> emittedCliqueKeys;

  vector<ui> eIndex;
  vector<ui> eIndexStamp;
  ui eIndexToken;

  void intersectInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  void setDiffInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  vector<ui> unionSet(const vector<ui> &A, const vector<ui> &B);

  vector<ui> commonExpand(const vector<ui> &E, const vector<ui> &S);
  vector<ui> collectAllCoveringCliques(const vector<ui> &M);
  vector<vector<ui>> buildHitSets(const vector<ui> &E,
                                  const vector<ui> &cliqueIds);
  vector<vector<ui>> singletonBranches(const vector<ui> &E);
  vector<vector<ui>>
  generateExactSiblingSets(const vector<ui> &E,
                           const vector<ui> &coveringCliqueIds,
                           bool *usePivotFallback = nullptr);
  vector<vector<ui>> efficientHittingSet(const vector<ui> &E,
                                         const vector<vector<ui>> &hitSets,
                                         bool *usePivotFallback = nullptr);
  void recordSolverCallStats(ui eSize, ui hSize);
  void recordSolverCompatStats(ull eligible, ull survivors);

  bool adj(ui u, ui v) const { return adjSet[u].count(v) != 0; }
  bool findOnePure(const vector<ui> &M, const vector<ui> &Q, vector<ui> &found);
  bool findOnePureRecursive(vector<ui> &R, vector<ui> P, vector<ui> X,
                            vector<ui> &found);
  ui pureNeighborsInP(ui u, const vector<ui> &P) const;
  void scanPurePXRState(const vector<ui> &P, const vector<ui> &X, ui &pivot,
                        ui &minPScore, ui &universalP, bool &xUniversal) const;
  void pureMatchingParts(const vector<ui> &P, vector<ui> &forced,
                         vector<pair<ui, ui>> &missingEdges) const;
  void enumerateAllPureBranch(const vector<ui> &M, const vector<ui> &Q);
  void enumerateAllPureBranchRecursive(vector<ui> &R, vector<ui> P,
                                       vector<ui> X);
  bool recordPureClique(vector<ui> C);

public:
  explicit ReorderSib(Graph &g, DegOrder order = DegOrder::ASCENDING,
                      ui minCliqueSize = 3);

  void findAllMaximalCliquesPure();
  void setSolverWorkBudget(ull budget) { solverWorkBudget = budget; }
  ull getCliqueCount() const { return cliqueCount; }
  ull getDuplicateCount() const { return dupBlocked; }
  ui getMaxCliqueSize() const { return static_cast<ui>(maxCliqueSize); }
  ull getChecksCount() const { return checksCount; }
  ull getBudgetFallbackCount() const { return solverBudgetFallbacks; }
  vector<vector<ui>> getCliques() const;
};
