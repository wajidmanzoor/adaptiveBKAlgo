#pragma once

#include "checked_count.h"
#include "common.h"
#include "graph.h"

enum class SibMethod {
  BACKTRACKING,
  OPTIMIZED
};

struct ReorderSibTestAccess;

// Optimized Adjacency List based Bron-Kerbosch with Pivoting and Pruning
class PivotBK {
private:
  ui n;
  vector<vector<ui>> adjList;
  ui cliqueCount;
  ui maxCliqueSize;
  ui checksCount;
  ui minCliqueSize;

  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);

public:
  explicit PivotBK(Graph &g, ui minCliqueSize = 3);

  void findAllMaximalCliques();
};

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
  vector<unordered_set<ui>> adjSet; // O(1) adjacency lookup
  ull cliqueCount;
  ull dupBlocked;
  size_t maxCliqueSize;
  ull externalCliqueCount;
  size_t externalMaxCliqueSize;
  ull checksCount;
  ull solverWorkBudget;
  ull solverBudgetFallbacks;
  SibMethod method;
  ui minCliqueSize;
  bool sp1;    // solver: unit propagation — force candidates that are the sole cover of a constraint
  bool sp2;    // solver: constraint subsumption — drop constraints implied by tighter ones
  bool sp3;    // solver: sort hitSets by ascending size so fail-first hits the hardest constraint first
  vector<vector<ui>> allCliques;
  // The search and cover index use internal permuted labels. External
  // validation restores original graph IDs through this inverse permutation.
  vector<ui> internalToOriginal;
  // Clique index: [vertex][level] = list of clique IDs. Pure mode stores its
  // emitted cliques in level zero and uses the index for cover lookup.
  vector<vector<vector<ui>>> cliquesByVertexByLevel;
  vector<ull> cliqueCountByVertex; // total cliques per vertex — for seed selection
  unordered_set<string> emittedCliqueKeys;

  vector<ui> eIndex;             // reusable vertex -> local E index map
  vector<ui> eIndexStamp;        // stamp for entries valid in current solver call
  ui eIndexToken;

  vector<ui> intersect(const vector<ui> &A, const vector<ui> &B);
  void intersectInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                              const vector<ui> &B,
                              const vector<ui> &exclude);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  void setDiffInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  vector<ui> unionSet(const vector<ui> &A, const vector<ui> &B);

  bool hitsAll(const vector<ui> &S, const vector<vector<ui>> &hitSets);
  vector<ui> commonExpand(const vector<ui> &E, const vector<ui> &S);
  vector<vector<ui>> buildHitSets(const vector<ui> &E,
                                  const vector<ui> &cliqueIds,
                                  ui maxHitSets = UINT_MAX);
  vector<vector<ui>> singletonBranches(const vector<ui> &E);

  vector<vector<ui>> backtrackingBranchBound(const vector<ui> &E,
                                             const vector<vector<ui>> &hitSets);
  vector<vector<ui>> efficientHittingSet(const vector<ui> &E,
                                         const vector<vector<ui>> &hitSets,
                                         bool *usePivotFallback = nullptr);
  void recordSolverCallStats(ui eSize, ui hSize);
  void recordSolverCompatStats(ull eligible, ull survivors);
  bool adj(ui u, ui v) const { return adjSet[u].count(v); }

  vector<ui> collectAllCoveringCliques(const vector<ui> &M);
  vector<vector<ui>>
  generateExactSiblingSets(const vector<ui> &E,
                           const vector<ui> &coveringCliqueIds,
                           bool *usePivotFallback = nullptr);
  bool findOnePure(const vector<ui> &M, const vector<ui> &Q,
                   vector<ui> &found);
  bool findOnePureRecursive(vector<ui> &R, vector<ui> P, vector<ui> X,
                            vector<ui> &found);
  ui pureNeighborsInP(ui u, const vector<ui> &P) const;
  void scanPurePXRState(const vector<ui> &P, const vector<ui> &X,
                        ui &pivot, ui &minPScore, ui &universalP,
                        bool &xUniversal) const;
  void pureMatchingParts(const vector<ui> &P, vector<ui> &forced,
                         vector<pair<ui, ui>> &missingEdges) const;
  void enumerateAllPureBranch(const vector<ui> &M, const vector<ui> &Q);
  void enumerateAllPureBranchRecursive(vector<ui> &R, vector<ui> P,
                                       vector<ui> X);
  bool recordPureClique(vector<ui> C);

public:
  ReorderSib(Graph &g, SibMethod method = SibMethod::OPTIMIZED,
             ui hitSetLimit = UINT_MAX,
             bool prune1 = true, bool prune2 = true,
             bool sp1 = true, bool sp2 = true, bool sp3 = true,
             bool sp4 = true, bool sp5 = true, bool sp6 = true,
             ui minCliqueSize = 3);
  void findAllMaximalCliquesPure();
  void setExternalResults(ull count, size_t maximumSize) {
    externalCliqueCount = count;
    externalMaxCliqueSize = maximumSize;
  }
  void setSolverWorkBudget(ull budget) { solverWorkBudget = budget; }
  ull getCliqueCount() const { return cliqueCount; }
  ull getDuplicateCount() const { return dupBlocked; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
  vector<vector<ui>> getCliques() const;
};
