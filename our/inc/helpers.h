#pragma once

#include "checked_count.h"
#include "common.h"
#include "graph.h"

enum class DegOrder { ORIGINAL, ASCENDING, DESCENDING };
enum class SibMethod {
  BACKTRACKING,
  OPTIMIZED
};

ui graphDegeneracy(const Graph &g);
struct ReorderSibTestAccess;

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

class BitsetBK {
private:
  ui n;
  ui words;
  ui degeneracy;
  bool lowDegreeGraph;
  ull lowDegreeCliqueCount;
  ui lowDegreeMaxSize;
  bool bipartite;
  bool chordalGraph;
  ull chordalCliqueCount;
  ui chordalMaxSize;
  bool twinModuleQuotient;
  ull twinModuleQuotientCount;
  ui twinModuleQuotientMaxSize;
  bool falseTwinQuotient;
  ull falseTwinQuotientCount;
  ui falseTwinQuotientMaxSize;
  bool trueTwinQuotient;
  ull trueTwinQuotientCount;
  ui trueTwinQuotientMaxSize;
  bool cliqueComponents;
  ull cliqueComponentCount;
  ui cliqueComponentMaxSize;
  bool completeMultipartite;
  ui completeMultipartiteParts;
  ull completeMultipartiteCount;
  vector<ull> adjBits;
  vector<ui> order;
  vector<ui> position;
  vector<vector<ull>> depthP;
  vector<vector<ull>> depthX;
  vector<vector<ull>> depthCand;
  vector<vector<ui>> depthActive;
  ull cliqueCount;
  ui maxCliqueSize;
  ull checksCount;

  const ull *neighbors(ui v) const;
  void ensureDepth(ui depth);
  bool isConnected(ui u, ui v) const;
  bool solveTwinModuleQuotient(ull &quotientCount, ui &quotientMaxSize) const;
  bool solveFalseTwinQuotient(ull &quotientCount, ui &quotientMaxSize) const;
  bool solveTrueTwinQuotient(ull &quotientCount, ui &quotientMaxSize) const;
  void detectCompleteMultipartite();
  bool isEmpty(const vector<ull> &bits, const vector<ui> &active) const;
  bool hasEdgeInP(const vector<ull> &P, const vector<ui> &active) const;
  bool tryComplementMatching(const vector<ull> &P, const vector<ui> &active,
                             ui rSize, ui pSize);
  bool tryComplementDegreeTwo(const vector<ull> &P, const vector<ull> &X,
                              const vector<ui> &active, ui rSize,
                              ui pSize);
  ui choosePivot(const vector<ull> &P, const vector<ull> &X,
                 const vector<ui> &active, ui &pSize, int &minPScore,
                 bool &xExtendsP) const;
  void bronKerboschRecursive(ui rSize, ui depth);

public:
  BitsetBK(Graph &g);
  void findAllMaximalCliques();
};

class LocalBitsetBK {
private:
  struct NeighborRange {
    const ui *first;
    const ui *last;

    const ui *begin() const { return first; }
    const ui *end() const { return last; }
    size_t size() const { return static_cast<size_t>(last - first); }
  };

  ui n;
  ui degeneracy;
  bool lowDegreeGraph;
  ull lowDegreeCliqueCount;
  ui lowDegreeMaxSize;
  bool bipartite;
  bool chordalGraph;
  ull chordalCliqueCount;
  ui chordalMaxSize;
  bool twinModuleQuotient;
  ull twinModuleQuotientCount;
  ui twinModuleQuotientMaxSize;
  bool cliqueComponents;
  ull cliqueComponentCount;
  ui cliqueComponentMaxSize;
  const Graph *graph;
  vector<ui> order;
  vector<ui> position;

  vector<ui> localVerts;
  ui localSize;
  ui localWords;
  vector<ull> localAdjBits;
  vector<ui> localIndex;
  vector<ui> localStamp;
  ui localToken;
  vector<ui> xStamp;
  ui xToken;

  vector<vector<ull>> depthP;
  vector<vector<ui>> depthX;
  vector<vector<ull>> depthCand;

  ull cliqueCount;
  ui maxCliqueSize;
  ull checksCount;

  void ensureDepth(ui depth);
  bool shouldDetectCliqueComponents(const Graph &g) const;
  bool detectCliqueComponents(const Graph &g);
  NeighborRange adj(ui u) const;
  bool isConnected(ui u, ui v) const;
  const ull *localNeighbors(ui idx) const;
  bool isEmpty(const vector<ull> &bits) const;
  ui popcount(const vector<ull> &bits) const;
  ui choosePivot(const vector<ull> &P, const vector<ui> &X, ui &pSize,
                 int &minPScore, bool &xExtendsP) const;
  bool xExtendsP(const vector<ui> &X, const vector<ull> &P,
                 ui pSize) const;
  void bronKerboschRecursive(ui rSize, ui depth);
  bool buildRoot(ui root);

public:
  LocalBitsetBK(Graph &g);
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
  // The search and cover index use internal permuted labels. External
  // validation restores original graph IDs through this inverse permutation.
  vector<ui> internalToOriginal;
  // Two-tier clique index: [vertex][level] = list of clique IDs.
  // Replaces the flat cliquesByVertex + foundLevel pair; with prune2 on we
  // jump directly to the relevant level bucket instead of scanning all cliques.
  vector<vector<vector<ui>>> cliquesByVertexByLevel;
  vector<ull> cliqueCountByVertex; // total cliques per vertex — for seed selection
  unordered_set<string> claimedEmptyBranches;
  unordered_set<string> emittedCliqueKeys;

  vector<char> lab; // lab[v] = 1 iff v is in current P; 2 iff in X

  // Pre-allocated per-depth workspace — avoids heap allocation inside enumerate.
  // Each depth slot holds the working P, X, and scratch P_new/X_new buffers.
  static constexpr ui MAX_ENUM_DEPTH = 200;
  ui enumDepth;
  vector<vector<ui>> depthPal;   // P_at_level per depth
  vector<vector<ui>> depthLX;    // localX per depth
  vector<vector<ui>> depthPnew;  // P_new scratch per depth
  vector<vector<ui>> depthXnew;  // X_new scratch per depth
  vector<ui> eIndex;             // reusable vertex -> local E index map
  vector<ui> eIndexStamp;        // stamp for entries valid in current solver call
  ui eIndexToken;
  vector<ui> collectVertexStamp; // reusable membership marks for collect inputs
  ui collectVertexToken;
  vector<ui> collectCliqueStamp; // reusable clique-ID marks for seed intersection
  ui collectCliqueToken;

  vector<ui> intersect(const vector<ui> &A, const vector<ui> &B);
  void intersectInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                              const vector<ui> &B,
                              const vector<ui> &exclude);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  void setDiffInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  vector<ui> unionSet(const vector<ui> &A, const vector<ui> &B);
  void unionInto(vector<ui> &out, const vector<ui> &A, const vector<ui> &B);
  void insertSortedVertex(vector<ui> &out, const vector<ui> &base, ui v);

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
                                         const vector<vector<ui>> &hitSets,
                                         bool *usePivotFallback = nullptr);
  void recordSolverCallStats(ui eSize, ui hSize);
  void recordSolverCompatStats(ull eligible, ull survivors);
  bool branchSpaceInsideClique(const vector<ui> &M, const vector<ui> &E,
                               const vector<ui> &C);
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

  void commitCliqueAndReorder(vector<ui> C,
                              vector<vector<ui>> &mustin,
                              vector<vector<ui>> &expandTo,
                              vector<char> &fullSkipCheck, ui treeIndex,
                              ui level, bool &done);

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
             bool sp4 = true, bool sp5 = true, bool sp6 = true,
             ui minCliqueSize = 3);
  void findAllMaximalCliques();
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
