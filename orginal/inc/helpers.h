#pragma once

#include "common.h"
#include "graph.h"

enum class DegOrder { ORIGINAL, ASCENDING, DESCENDING };
enum class SibMethod {
  BACKTRACKING,
  OPTIMIZED
};

// Optimized Adjacency List based Bron-Kerbosch with Pivoting and Pruning
class PivotBK {
private:
  ui n;
  vector<vector<ui>> adjList;
  ui cliqueCount;
  ui maxCliqueSize;
  ui checksCount;

  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);

public:
  PivotBK(Graph &g, DegOrder order = DegOrder::ASCENDING);

  void findAllMaximalCliques();
};

class ReorderSib {
private:
  ui n;
  vector<vector<ui>> adjList;
  vector<vector<ui>> adjList2;
  vector<unordered_set<ui>> adjSet; // O(1) adjacency lookup
  ui cliqueCount;
  ui dupBlocked;
  size_t maxCliqueSize;
  ui checksCount;
  SibMethod method;
  ui hitSetLimit;
  bool prune1; // dominance pruning: drop Ci when Ci ⊂ Cj among covering cliques
  bool prune2; // level filter: skip covering cliques older than level-1
  // prune3 is always ON: aggressive branch skip via branchSpaceInsideClique + fullSkipCheck
  bool sp1;    // solver: unit propagation — force candidates that are the sole cover of a constraint
  bool sp2;    // solver: constraint subsumption — drop constraints implied by tighter ones
  bool sp3;    // solver: sort hitSets by ascending size so fail-first hits the hardest constraint first
  bool sp4;    // collect: skip covering cliques whose intersection with E is exactly M (trivially weak)
  bool sp5;    // rCall: skip branches where |mustin|+|expandTo| <= 2 (can't form clique of size > 2)
  bool sp6;    // enumerate: skip reordered branches with empty expandTo and |mustin| <= 2
  vector<vector<ui>> allCliques;
  // Two-tier clique index: [vertex][level] = list of clique IDs.
  // Replaces the flat cliquesByVertex + foundLevel pair; with prune2 on we
  // jump directly to the relevant level bucket instead of scanning all cliques.
  vector<vector<vector<ui>>> cliquesByVertexByLevel;
  vector<ui> cliqueCountByVertex; // total cliques per vertex — for seed selection
  unordered_set<string> claimedEmptyBranches;

  vector<char> lab; // lab[v] = 1 iff v is in current P; 2 iff in X

  // Pre-allocated per-depth workspace — avoids heap allocation inside enumerate.
  // Each depth slot holds the working P, X, and scratch P_new/X_new buffers.
  static constexpr ui MAX_ENUM_DEPTH = 200;
  ui enumDepth;
  vector<vector<ui>> depthPal;   // P_at_level per depth
  vector<vector<ui>> depthLX;    // localX per depth
  vector<vector<ui>> depthPnew;  // P_new scratch per depth
  vector<vector<ui>> depthXnew;  // X_new scratch per depth

  vector<ui> intersect(const vector<ui> &A, const vector<ui> &B);
  void intersectInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                              const vector<ui> &B,
                              const vector<ui> &exclude);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  void setDiffInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  vector<ui> unionSet(const vector<ui> &A, const vector<ui> &B);
  void unionInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);

  bool hitsAll(const vector<ui> &S, const vector<vector<ui>> &hitSets);
  vector<ui> commonExpand(const vector<ui> &E, const vector<ui> &S);
  vector<ui> collectCoveringCliques(const vector<ui> &M, const vector<ui> &E, ui level);
  vector<vector<ui>> buildHitSets(const vector<ui> &E,
                                  const vector<ui> &cliqueIds,
                                  ui maxHitSets = UINT_MAX);
  vector<vector<ui>> singletonBranches(const vector<ui> &E);
  vector<vector<ui>> minimalByInclusion(vector<vector<ui>> solutions);

  vector<ui> pruneByDominance(const vector<ui> &cliqueIds);
  vector<vector<ui>>
  generateSiblingSetsFromCliques(const vector<ui> &E,
                                 const vector<ui> &cliqueIds);
  vector<vector<ui>> backtrackingBranchBound(const vector<ui> &E,
                                             const vector<vector<ui>> &hitSets);
  vector<vector<ui>> efficientHittingSet(const vector<ui> &E,
                                         const vector<vector<ui>> &hitSets);
  bool branchSpaceInsideClique(const vector<ui> &M, const vector<ui> &E,
                               const vector<ui> &C);
  bool adj(ui u, ui v) const { return adjSet[u].count(v); }

  void rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo, ui level,
             vector<char> fullSkipCheck);
  void enumerate(vector<ui> &R, const vector<ui> &P, const vector<ui> &X,
                 vector<vector<ui>> &mustin, vector<vector<ui>> &expandTo,
                 vector<char> &fullSkipCheck, ui treeIndex, ui level,
                 bool &done);

public:
  ReorderSib(Graph &g, DegOrder order = DegOrder::ORIGINAL,
             SibMethod method = SibMethod::OPTIMIZED,
             ui hitSetLimit = UINT_MAX,
             bool prune1 = true, bool prune2 = true,
             bool sp1 = true, bool sp2 = true, bool sp3 = true,
             bool sp4 = true, bool sp5 = true, bool sp6 = true);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};
