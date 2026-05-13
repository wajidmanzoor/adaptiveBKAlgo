#pragma once

#include "common.h"
#include "graph.h"

enum class DegOrder { ORIGINAL, ASCENDING, DESCENDING };
enum class SibMethod {
  BRUTE_FORCE,
  BACKTRACKING,
  GREEDY,
  BITMASK,
  MIN_HITTING_SET,
  OPTIMIZED
};

// Optimized Adjacency List based Bron-Kerbosch with Pivoting and Pruning
class PivotBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  ui maxCliqueSize; // for neighborhood size pruning
  ui checksCount;   // total vertex-set checks in bronKerboschRecursive
  ui redundancy;

  vector<vector<ui>> redendantChecks; // Track redundant vertices for pruning

  vector<vector<ui>>
      foundCliques; // Store found maximal cliques for output/debugging

  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);
  bool isPSubsetOfFoundClique(const vector<ui> &P);

public:
  PivotBK(Graph &g, DegOrder order = DegOrder::ASCENDING);

  void findAllMaximalCliques();
};

class ReorderSib {
private:
  ui n;
  vector<vector<ui>> adjList;
  vector<vector<ui>> adjList2;
  ui cliqueCount;
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
  vector<ui> foundLevel;
  vector<vector<ui>> cliquesByVertex;
  unordered_set<string> seenCliques;

  vector<ui> intersect(const vector<ui> &A, const vector<ui> &B);
  vector<ui> setDiff(const vector<ui> &A, const vector<ui> &B);
  vector<ui> unionSet(const vector<ui> &A, const vector<ui> &B);
  vector<ui> compliment(const vector<ui> &vector1);

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
  vector<vector<ui>> bruteForceBySize(const vector<ui> &E,
                                      const vector<vector<ui>> &hitSets);
  vector<vector<ui>> backtrackingBranchBound(const vector<ui> &E,
                                             const vector<vector<ui>> &hitSets);
  vector<vector<ui>> greedyApproximation(const vector<ui> &E,
                                         const vector<vector<ui>> &hitSets);
  vector<vector<ui>> bitmaskExactSearch(const vector<ui> &E,
                                        const vector<vector<ui>> &hitSets);
  vector<vector<ui>> minimumCliqueHittingSet(const vector<ui> &E,
                                             const vector<vector<ui>> &hitSets);
  vector<vector<ui>> efficientHittingSet(const vector<ui> &E,
                                         const vector<vector<ui>> &hitSets);
  bool branchSpaceInsideClique(const vector<ui> &M, const vector<ui> &E,
                               const vector<ui> &C);

  void rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo, ui level,
             vector<char> fullSkipCheck);
  void enumerate(vector<ui> &R, vector<ui> &Q, vector<vector<ui>> &mustin,
                 vector<vector<ui>> &expandTo, vector<char> &fullSkipCheck,
                 ui treeIndex, ui level, bool &done);

public:
  ReorderSib(Graph &g, DegOrder order = DegOrder::ORIGINAL,
             SibMethod method = SibMethod::BRUTE_FORCE,
             ui hitSetLimit = UINT_MAX,
             bool prune1 = true, bool prune2 = true,
             bool sp1 = true, bool sp2 = true, bool sp3 = true,
             bool sp4 = true, bool sp5 = true, bool sp6 = true);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};
