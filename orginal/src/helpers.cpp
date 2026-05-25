#include "../inc/helpers.h"
#include <chrono>
#include <functional>
#include <iomanip>
#include <numeric>
// Returns peelSeq index of verticies by core value
// peelSeq[0] = highest-core vertex, peelSeq[n-1] = lowest.
static vector<ui> computePeelSeq(const Graph &g) {
  ui n = g.n;
  vector<ui> deg(g.degree.begin(), g.degree.end());
  ui maxDeg = *max_element(deg.begin(), deg.end());

  vector<ui> bins(maxDeg + 1, 0);
  for (ui d : deg)
    bins[d]++;
  vector<ui> binStart(maxDeg + 1, 0);
  partial_sum(bins.begin(), bins.end() - 1, binStart.begin() + 1);

  vector<ui> pos(n), sorted(n);
  for (ui v = 0; v < n; v++) {
    pos[v] = binStart[deg[v]]++;
    sorted[pos[v]] = v;
  }
  for (ui d = 0; d <= maxDeg; d++)
    binStart[d] -= bins[d]; // reset to start of each bin

  vector<ui> peelSeq(n);
  for (ui i = 0; i < n; i++) {
    ui v = sorted[i];
    peelSeq[n - 1 - i] = v; // ascending peel order → descending in array
    for (ui j = g.offset[v]; j < g.offset[v + 1]; j++) {
      ui u = g.neighbors[j];
      if (deg[u] > deg[v]) {
        ui du = deg[u];
        ui pu = pos[u];
        ui pw = binStart[du];
        ui w = sorted[pw];
        if (u != w) {
          pos[u] = pw;
          sorted[pu] = w;
          pos[w] = pu;
          sorted[pw] = u;
        }
        binStart[du]++;
        deg[u]--;
      }
    }
  }
  return peelSeq;
}

// Build adjList (all neighbors) and adjList2 (only neighbors with higher index
// in the permuted order).
static void buildAdjLists(const Graph &g, const vector<ui> &perm,
                          vector<vector<ui>> &adjList,
                          vector<vector<ui>> &adjList2) {
  ui n = g.n;
  adjList.assign(n, {});
  adjList2.assign(n, {});
  for (ui u = 0; u < n; u++) {
    ui nu = perm[u];
    for (ui j = g.offset[u]; j < g.offset[u + 1]; j++) {
      ui nv = perm[g.neighbors[j]];
      adjList[nu].push_back(nv);
      if (nv > nu)
        adjList2[nu].push_back(nv);
    }
    sort(adjList[nu].begin(), adjList[nu].end());
    sort(adjList2[nu].begin(), adjList2[nu].end());
  }
}

// PivotBK Implementation
PivotBK::PivotBK(Graph &g, DegOrder order) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  vector<ui> perm(n);

  // Canonical order
  if (order == DegOrder::ORIGINAL) {
    for (ui i = 0; i < n; i++)
      perm[i] = i;
  } else {
    vector<ui> peelSeq = computePeelSeq(g);

    // Ascending degeneracy order
    if (order == DegOrder::ASCENDING) {
      for (ui i = 0; i < n; i++)
        perm[peelSeq[n - 1 - i]] = i; // low-core → low index
    } else {
      // Descending degeneracy order
      for (ui i = 0; i < n; i++)
        perm[peelSeq[i]] = i; // high-core → low index
    }
  }

  adjList.resize(n);
  vector<vector<ui>> dummy;
  buildAdjLists(g, perm, adjList, dummy);
}

vector<ui> PivotBK::intersect(const vector<ui> &set1,
                              const vector<ui> &neighbors) {
  vector<ui> result;
  result.reserve(min(set1.size(), neighbors.size()));
  ui i = 0, j = 0;
  while (i < set1.size() && j < neighbors.size()) {
    if (set1[i] == neighbors[j]) {
      result.push_back(set1[i]);
      i++;
      j++;
    } else if (set1[i] < neighbors[j]) {
      i++;
    } else {
      j++;
    }
  }
  return result;
}
bool PivotBK::isEmpty(const vector<ui> &set) { return set.empty(); }

bool PivotBK::isConnected(ui u, ui v) {
  // Binary search in sorted adjacency list
  // if v and u are connected.
  return binary_search(adjList[u].begin(), adjList[u].end(), v);
}

ui PivotBK::choosePivot(const vector<ui> &P, const vector<ui> &X) {
  ui bestPivot = P.empty() ? (X.empty() ? 0 : X[0]) : P[0];
  ui maxElimination = 0;

  // Check vertices in P to find the one that maximizes |P ∩ N(u)|
  for (ui u : P) {
    ui elimination = (ui)intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }

  // check vertices in X as well, since pivot can be from P ∪ X
  for (ui u : X) {
    ui elimination = (ui)intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }
  return bestPivot;
}
void PivotBK::bronKerboschRecursive(vector<ui> &R, vector<ui> &P,
                                    vector<ui> &X) {
  checksCount++;

  // Basic pruning: check if P and X are empty
  // Clique found
  if (isEmpty(P) && isEmpty(X)) {
    if (R.size() <= 2)
      return;

    // Found a maximal clique
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());
    return;
  }

  // P empty but X non-empty: R is NOT maximal, prune
  if (isEmpty(P))
    return;

  // Choose pivot from P ∪ X such that |P ∩ N(pivot)| is maximized to minimize
  // recursive calls
  ui pivot = choosePivot(P, X);

  // P = P \ N(pivot)
  vector<ui> candidates;
  candidates.reserve(P.size());
  for (ui v : P) {
    if (!isConnected(v, pivot))
      candidates.push_back(v);
  }

  // For each candidate vertex v, we add it to the growing clique R and
  // recursively explore.
  for (ui v : candidates) {
    // add v to partial clique R
    R.push_back(v);

    // new P =  P ∩ N(v), new X = X ∩ N(v)
    vector<ui> new_P = intersect(P, adjList[v]);
    vector<ui> new_X = intersect(X, adjList[v]);

    // recurse with new sets
    bronKerboschRecursive(R, new_P, new_X);

    R.pop_back(); // backtrack

    // Move v from P to X, keeping X sorted
    auto it = find(P.begin(), P.end(), v);
    if (it != P.end())
      P.erase(it);
    X.insert(lower_bound(X.begin(), X.end(), v), v);
  }
}
void PivotBK::findAllMaximalCliques() {
  vector<ui> R;
  vector<ui> X;
  vector<ui> P(n);

  for (ui i = 0; i < n; i++)
    P[i] = i;

  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  auto t0 = chrono::high_resolution_clock::now();
  bronKerboschRecursive(R, P, X);
  auto t1 = chrono::high_resolution_clock::now();
  double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
  cout << "Total Vertex-Set Checks: " << checksCount << endl;
  cout << fixed << setprecision(3) << "Time: " << ms << " ms" << endl;
}

// ReorderSib profiling
struct RSibProf {
  double rCall_ms = 0, enumerate_ms = 0, collect_ms = 0, solver_ms = 0;
  double minimal_ms = 0, commonExp_ms = 0, buildHit_ms = 0;
  long rCall_n = 0, enumerate_n = 0, collect_n = 0, solver_n = 0;
  long minimal_n = 0, commonExp_n = 0, buildHit_n = 0;
  void reset() { *this = RSibProf{}; }
  void print(double total_ms) const {
    auto pct = [&](double v) {
      return total_ms > 0 ? 100.0 * v / total_ms : 0.0;
    };
    printf(
        "\n── ReorderSib cost breakdown ─────────────────────────────────\n");
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "enumerate (core search)",
           enumerate_ms, pct(enumerate_ms), enumerate_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "rCall overhead", rCall_ms,
           pct(rCall_ms), rCall_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "solver", solver_ms,
           pct(solver_ms), solver_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "collectCoveringCliques",
           collect_ms, pct(collect_ms), collect_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "buildHitSets",
           buildHit_ms, pct(buildHit_ms), buildHit_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "minimalByInclusion",
           minimal_ms, pct(minimal_ms), minimal_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "commonExpand",
           commonExp_ms, pct(commonExp_ms), commonExp_n);
    printf("  %-30s %9.3f ms  %5.1f%%\n", "TOTAL (wall)", total_ms, 100.0);
    printf("──────────────────────────────────────────────────────────────\n");
  }
};
static RSibProf rsp;

#if PROFILING
struct ScopedTimer {
  chrono::high_resolution_clock::time_point t0;
  double &acc;
  long &cnt;
  ScopedTimer(double &a, long &c)
      : t0(chrono::high_resolution_clock::now()), acc(a), cnt(c) {}
  ~ScopedTimer() {
    acc += chrono::duration<double, milli>(
               chrono::high_resolution_clock::now() - t0)
               .count();
    ++cnt;
  }
};
#else
struct ScopedTimer {
  ScopedTimer(double &, long &) {}
};
#endif

// ReorderSib Implementation
ReorderSib::ReorderSib(Graph &g, DegOrder order, SibMethod method,
                       ui hitSetLimit, bool prune1, bool prune2, bool sp1,
                       bool sp2, bool sp3, bool sp4, bool sp5, bool sp6)
    : hitSetLimit(hitSetLimit), prune1(prune1), prune2(prune2), sp1(sp1),
      sp2(sp2), sp3(sp3), sp4(sp4), sp5(sp5), sp6(sp6) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  this->method = method;
  cliquesByVertexByLevel.resize(n);
  cliqueCountByVertex.resize(n, 0);

  vector<ui> perm(n);

  // Canonical order
  if (order == DegOrder::ORIGINAL) {
    for (ui i = 0; i < n; i++)
      perm[i] = i;
  } else {
    // Degeneracy order Ascending or Descending
    vector<ui> peelSeq = computePeelSeq(g);
    if (order == DegOrder::ASCENDING) {
      for (ui i = 0; i < n; i++)
        perm[peelSeq[n - 1 - i]] = i;
    } else {
      for (ui i = 0; i < n; i++)
        perm[peelSeq[i]] = i;
    }
  }

  buildAdjLists(g, perm, adjList, adjList2);

  adjSet.resize(n);
  for (ui u = 0; u < n; u++)
    for (ui v : adjList[u])
      adjSet[u].insert(v);

  lab.assign(n, 0);
  enumDepth = 0;
  const ui depthCap = max<ui>(MAX_ENUM_DEPTH, n + 1);
  depthPal.resize(depthCap);
  depthLX.resize(depthCap);
  depthPnew.resize(depthCap);
  depthXnew.resize(depthCap);
}

vector<ui> ReorderSib::intersect(const vector<ui> &A, const vector<ui> &B) {
  vector<ui> C;
  C.reserve(min(A.size(), B.size()));
  ui i = 0, j = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] == B[j]) {
      C.push_back(A[i]);
      i++;
      j++;
    } else if (A[i] < B[j])
      i++;
    else
      j++;
  }
  return C;
}

void ReorderSib::intersectInto(vector<ui> &out, const vector<ui> &A,
                               const vector<ui> &B) {
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  ui i = 0, j = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] == B[j]) {
      out.push_back(A[i]);
      i++;
      j++;
    } else if (A[i] < B[j]) {
      i++;
    } else {
      j++;
    }
  }
}

vector<ui> ReorderSib::setDiff(const vector<ui> &A, const vector<ui> &B) {
  vector<ui> C;
  C.reserve(A.size());
  ui i = 0, j = 0;
  while (i < A.size()) {
    if (j == (ui)B.size() || A[i] < B[j]) {
      C.push_back(A[i]);
      i++;
    } else if (A[i] == B[j]) {
      i++;
      j++;
    } else
      j++;
  }
  return C;
}

vector<ui> ReorderSib::unionSet(const vector<ui> &A, const vector<ui> &B) {
  vector<ui> U;
  U.reserve(A.size() + B.size());
  ui i = 0, j = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] < B[j]) {
      U.push_back(A[i]);
      i++;
    } else if (A[i] > B[j]) {
      U.push_back(B[j]);
      j++;
    } else {
      U.push_back(A[i]);
      i++;
      j++;
    }
  }
  while (i < A.size())
    U.push_back(A[i++]);
  while (j < B.size())
    U.push_back(B[j++]);
  return U;
}

bool ReorderSib::hitsAll(const vector<ui> &S,
                         const vector<vector<ui>> &hitSets) {
  for (const vector<ui> &hitSet : hitSets) {
    bool hit = false;
    for (ui v : S) {
      if (binary_search(hitSet.begin(), hitSet.end(), v)) {
        hit = true;
        break;
      }
    }
    if (!hit)
      return false;
  }
  return true;
}

// After choosing a sibling set S, only vertices still in E and adjacent to
// every vertex of S can continue to grow the branch.
vector<ui> ReorderSib::commonExpand(const vector<ui> &E, const vector<ui> &S) {
  ScopedTimer _t(rsp.commonExp_ms, rsp.commonExp_n);
  vector<ui> result = setDiff(E, S);
  for (ui v : S) {
    const auto &adj = adjList[v];
    ui k = 0, j = 0;
    for (ui i = 0; i < (ui)result.size(); i++) {
      while (j < (ui)adj.size() && adj[j] < result[i]) j++;
      if (j < (ui)adj.size() && adj[j] == result[i])
        result[k++] = result[i];
    }
    result.resize(k);
  }
  return result;
}

// Find previously discovered cliques that already contain M.
vector<ui> ReorderSib::collectCoveringCliques(const vector<ui> &M,
                                              const vector<ui> &E, ui level) {
  ScopedTimer _t(rsp.collect_ms, rsp.collect_n);
  vector<ui> result;
  if (M.empty())
    return result;

  // Seed: vertex in M whose total clique count is smallest (tightest filter).
  ui seed = M[0];
  for (ui v : M)
    if (cliqueCountByVertex[v] < cliqueCountByVertex[seed]) seed = v;

  // With prune2: only look at levels >= level-1 in the two-tier index,
  // skipping all older cliques without scanning them.
  const ui startLevel = prune2 ? (level > 0 ? level - 1 : 0) : 0;
  const auto &seedBuckets = cliquesByVertexByLevel[seed];
  for (ui l = startLevel; l < (ui)seedBuckets.size(); l++) {
    for (ui cId : seedBuckets[l]) {
      bool containsAllMustin = includes(
          allCliques[cId].begin(), allCliques[cId].end(), M.begin(), M.end());
      if (!containsAllMustin)
        continue;

      // sp4: skip cliques whose intersection with E is exactly M.
      if (sp4) {
        const auto &C = allCliques[cId];
        bool hasVertexBeyondM = false;
        for (ui v : C) {
          if (!binary_search(M.begin(), M.end(), v) &&
              binary_search(E.begin(), E.end(), v)) {
            hasVertexBeyondM = true;
            break;
          }
        }
        if (!hasVertexBeyondM)
          continue;
      }

      result.push_back(cId);
    }
  }
  return result;
}

// Each covering clique C contributes the constraint "pick something from E
// that is outside C", which is exactly E \ C.
vector<vector<ui>> ReorderSib::buildHitSets(const vector<ui> &E,
                                            const vector<ui> &cliqueIds,
                                            ui maxHitSets) {
  ScopedTimer _t(rsp.buildHit_ms, rsp.buildHit_n);

  // When capping, keep the cliques with the most overlap with E — those produce
  // the smallest hit sets (E \ C), which are the tightest constraints.
  // Tighter constraints → solver has fewer candidates per constraint → runs
  // faster.
  const vector<ui> *ids = &cliqueIds;
  vector<ui> sorted;
  if (cliqueIds.size() > maxHitSets) {
    sorted = cliqueIds;
    // Sort by clique size descending: larger clique → smaller E\C → tighter
    // constraint
    sort(sorted.begin(), sorted.end(), [&](ui a, ui b) {
      return allCliques[a].size() > allCliques[b].size();
    });
    sorted.resize(maxHitSets);
    ids = &sorted;
  }

  vector<vector<ui>> hitSets;
  hitSets.reserve(ids->size());
  for (ui cId : *ids)
    hitSets.push_back(setDiff(E, allCliques[cId]));
  return hitSets;
}

// If M is not covered by any old clique, the sibling effect does nothing and
// we branch on one vertex at a time exactly like the base reorder search.
vector<vector<ui>> ReorderSib::singletonBranches(const vector<ui> &E) {
  vector<vector<ui>> branches;
  branches.reserve(E.size());
  for (ui v : E)
    branches.push_back({v});
  return branches;
}

// Remove duplicate sibling sets, then keep only the ones that are minimal by
// inclusion. If S1 contains S2, S1 is unnecessary because S2 already enforces
// the same sibling split with a smaller branch seed.
vector<vector<ui>>
ReorderSib::minimalByInclusion(vector<vector<ui>> solutions) {
  ScopedTimer _t(rsp.minimal_ms, rsp.minimal_n);
  for (vector<ui> &S : solutions)
    sort(S.begin(), S.end());
  sort(solutions.begin(), solutions.end());
  solutions.erase(unique(solutions.begin(), solutions.end()), solutions.end());

  vector<vector<ui>> minimal;
  for (ui i = 0; i < solutions.size(); i++) {
    bool hasSmallerCover = false;
    for (ui j = 0; j < solutions.size(); j++) {
      if (i == j || solutions[j].size() >= solutions[i].size())
        continue;
      if (includes(solutions[i].begin(), solutions[i].end(),
                   solutions[j].begin(), solutions[j].end())) {
        hasSmallerCover = true;
        break;
      }
    }
    if (!hasSmallerCover)
      minimal.push_back(solutions[i]);
  }
  return minimal;
}

// Drop covering cliques that are proper subsets of another covering clique.
// If Ci ⊂ Cj then E\Ci ⊃ E\Cj, so Cj's hit-set constraint implies Ci's.
// Any solution that satisfies Cj automatically satisfies Ci, so Ci is
// redundant. Cliques in allCliques are stored sorted, so std::includes works
// directly.
vector<ui> ReorderSib::pruneByDominance(const vector<ui> &cliqueIds) {
  vector<ui> result;
  const ui k = (ui)cliqueIds.size();
  for (ui i = 0; i < k; i++) {
    const auto &Ci = allCliques[cliqueIds[i]];
    bool dominated = false;
    for (ui j = 0; j < k && !dominated; j++) {
      if (i == j)
        continue;
      const auto &Cj = allCliques[cliqueIds[j]];
      if (Cj.size() > Ci.size() &&
          includes(Cj.begin(), Cj.end(), Ci.begin(), Ci.end()))
        dominated = true;
    }
    if (!dominated)
      result.push_back(cliqueIds[i]);
  }
  return result;
}

// Convert the sibling effect into a clique-constrained hitting-set problem
// over the already discovered covering cliques, then solve it with the
// requested strategy.
vector<vector<ui>>
ReorderSib::generateSiblingSetsFromCliques(const vector<ui> &E,
                                           const vector<ui> &cliqueIds) {
  if (cliqueIds.empty())
    return singletonBranches(E);

  vector<ui> prunedVec;
  if (prune1)
    prunedVec = pruneByDominance(cliqueIds);
  const vector<ui> &pruned = prune1 ? prunedVec : cliqueIds;

  vector<vector<ui>> hitSets = buildHitSets(E, pruned, hitSetLimit);

  // An empty hit set means some old clique fully contains E, so no sibling
  // choice can separate the current branch from that clique.
  for (const vector<ui> &hitSet : hitSets) {
    if (hitSet.empty())
      return {};
  }

  if (pruned.size() == 1)
    return singletonBranches(hitSets[0]);

  if (method == SibMethod::BACKTRACKING)
    return backtrackingBranchBound(E, hitSets);

  return efficientHittingSet(E, hitSets);
}

vector<vector<ui>>
ReorderSib::backtrackingBranchBound(const vector<ui> &E,
                                    const vector<vector<ui>> &hitSets) {
  ScopedTimer _t(rsp.solver_ms, rsp.solver_n);
  vector<vector<ui>> solutions;
  vector<ui> current;

  // DFS over clique-compatible subsets of E. Once the current set already hits
  // every constraint, record it and stop descending that branch.
  function<void(ui)> dfs = [&](ui start) {
    if (hitsAll(current, hitSets)) {
      solutions.push_back(current);
      return;
    }

    for (ui i = start; i < E.size(); i++) {
      bool connected = true;
      for (ui v : current) {
        if (!binary_search(adjList[v].begin(), adjList[v].end(), E[i])) {
          connected = false;
          break;
        }
      }
      if (!connected)
        continue;

      current.push_back(E[i]);
      dfs(i + 1);
      current.pop_back();
    }
  };

  dfs(0);
  return minimalByInclusion(solutions);
}

// Optimized exact solver for all minimal clique-constrained hitting sets.
//
// Improvements over backtrackingBranchBound:
//   1. Bitmask coverage  — done-check and update are O(1) bitwise ops.
//   2. Incremental candidate list — compat-filtered frontier passed down,
//      no per-step binary_search into adjList.
//   3. Fail-first dead-branch  — prune as soon as any uncovered constraint
//      has zero candidates left.
//   4. "Covers nothing new" skip — a vertex that adds no new coverage can
//      never be part of a minimal solution; skip it unconditionally.
//   5. Live minimal-set maintenance — dominated solutions are removed the
//      moment a smaller one is found; no post-pass minimalByInclusion needed.
//   6. Coverage-descending candidate order — high-utility vertices tried
//      first, producing solutions earlier and enabling more pruning.
vector<vector<ui>>
ReorderSib::efficientHittingSet(const vector<ui> &E,
                                const vector<vector<ui>> &hitSets) {
  ScopedTimer _t(rsp.solver_ms, rsp.solver_n);
  const ui eSize = (ui)E.size();
  const ui hSize = (ui)hitSets.size();

  // Fall back if bitmask cannot cover all constraints.
  if (hSize > 63)
    return backtrackingBranchBound(E, hitSets);

  // ── Pre-computation ───────────────────────────────────────────────────────

  // compat[i*eSize+j] = 1  iff  E[i] and E[j] are adjacent in the graph.
  vector<char> compat(eSize * eSize, 0);
  for (ui i = 0; i < eSize; i++)
    for (ui j = i + 1; j < eSize; j++)
      if (binary_search(adjList[E[i]].begin(), adjList[E[i]].end(), E[j]))
        compat[i * eSize + j] = compat[j * eSize + i] = 1;

  // cov[i] = bitmask of hitSets that E[i] covers.
  // Non-const so sp2 can reduce it by clearing subsumed constraint bits.
  ull fullMask = (hSize == 64) ? ~0ULL : ((1ULL << hSize) - 1);
  vector<ull> cov(eSize, 0);
  {
    // ── sp3: Sort hitSets ascending by size before encoding into bitmasks. ──
    // Smaller hitSets are harder to satisfy (fewer candidates), so assigning
    // them to lower bit positions makes the fail-first dead-branch check
    // encounter the hardest constraint first and prune sooner.
    vector<ui> hOrder(hSize);
    iota(hOrder.begin(), hOrder.end(), 0);
    if (sp3)
      sort(hOrder.begin(), hOrder.end(),
           [&](ui a, ui b) { return hitSets[a].size() < hitSets[b].size(); });

    vector<ui> eIdx(n, eSize);
    for (ui i = 0; i < eSize; i++)
      eIdx[E[i]] = i;
    for (ui bit = 0; bit < hSize; bit++) {
      ui h = hOrder[bit];
      for (ui v : hitSets[h])
        if (eIdx[v] < eSize)
          cov[eIdx[v]] |= (1ULL << bit);
    }
  }

  // Initial candidate order: descending coverage count so high-utility
  // vertices are tried first, finding solutions sooner for better pruning.
  vector<ui> initCands(eSize);
  iota(initCands.begin(), initCands.end(), 0);
  sort(initCands.begin(), initCands.end(), [&](ui a, ui b) {
    return __builtin_popcountll(cov[a]) > __builtin_popcountll(cov[b]);
  });

  // ── sp1: Unit Propagation ─────────────────────────────────────────────────
  // For any constraint covered by exactly one candidate, that candidate is
  // forced into every solution. Pre-select all forced candidates, filter the
  // remaining candidates to be clique-compatible with them, and start the DFS
  // with the forced coverage already accumulated.
  vector<ui> forcedIdxs; // E-indices forced into every solution
  ull forcedCov = 0;

  if (sp1) {
    vector<ui> activeCands = initCands;
    bool changed = true;
    bool conflict = false;

    while (changed && !conflict) {
      changed = false;
      ull uncov = fullMask & ~forcedCov;
      if (!uncov)
        break;

      ull tmp = uncov;
      while (tmp && !conflict) {
        int h = __builtin_ctzll(tmp);
        tmp &= tmp - 1;

        ui sole = eSize;
        int cnt = 0;
        for (ui ci : activeCands) {
          if (cov[ci] & (1ULL << h)) {
            sole = ci;
            if (++cnt > 1)
              break;
          }
        }

        if (cnt == 0) {
          conflict = true;
          break;
        }
        if (cnt == 1) {
          // Check clique compatibility with already-forced vertices.
          for (ui fv : forcedIdxs) {
            if (!compat[fv * eSize + sole]) {
              conflict = true;
              break;
            }
          }
          if (conflict)
            break;

          forcedIdxs.push_back(sole);
          forcedCov |= cov[sole];

          vector<ui> next;
          next.reserve(activeCands.size());
          for (ui ci : activeCands)
            if (ci != sole && compat[sole * eSize + ci])
              next.push_back(ci);
          activeCands = std::move(next);
          changed = true;
          break; // restart scan with updated candidates
        }
      }
    }

    if (conflict)
      return {};

    initCands = std::move(activeCands);
  }

  // ── sp2: Constraint Subsumption ───────────────────────────────────────────
  // Constraint i is subsumed by constraint j when every candidate covering j
  // also covers i (coverSet(j) ⊆ coverSet(i)). Satisfying j then implies
  // satisfying i, so i can be dropped from fullMask.
  // We work on the post-sp1 initCands so forced-vertex coverage is reflected.
  if (sp2) {
    for (ui h = 0; h < hSize; h++) {
      if (!(fullMask & (1ULL << h)))
        continue; // already dropped
      if (forcedCov & (1ULL << h))
        continue; // h already covered by forced, DFS won't see it
      for (ui g = 0; g < hSize; g++) {
        if (g == h || !(fullMask & (1ULL << g)))
          continue;
        // g must not be force-covered: if no initCand covers g (because a
        // forced vertex was the sole cover and sp1 removed it), the subsumption
        // check would pass vacuously and incorrectly drop h.
        if (forcedCov & (1ULL << g))
          continue;
        // Check coverSet(g, initCands) ⊆ coverSet(h, initCands):
        // no remaining candidate covers g but not h.
        bool gSubsumesH = true;
        for (ui ci : initCands) {
          if ((cov[ci] & (1ULL << g)) && !(cov[ci] & (1ULL << h))) {
            gSubsumesH = false;
            break;
          }
        }
        if (gSubsumesH) {
          fullMask &= ~(1ULL << h); // drop h — implied by g
          break;
        }
      }
    }
  }

  // Early exit: if forced coverage (sp1) already satisfies all remaining
  // constraints (possibly reduced by sp2), return the forced solution directly.
  if ((forcedCov & fullMask) == fullMask) {
    if (forcedIdxs.empty())
      return {}; // nothing forced, nothing to cover
    vector<ui> sol;
    for (ui idx : forcedIdxs)
      sol.push_back(E[idx]);
    sort(sol.begin(), sol.end());
    return {sol};
  }

  // ── DFS ───────────────────────────────────────────────────────────────────

  // solutions is maintained as a live minimal-by-inclusion set throughout.
  vector<vector<ui>> solutions;
  // cur holds E-indices of the partial solution, in strictly increasing order.
  vector<ui> cur;

  // cands : E-indices still reachable (clique-compatible with cur, index >
  //         last element of cur). Passed by value so each level owns its copy.
  // covered: bitmask of constraints already satisfied by cur.
  function<void(vector<ui>, ull)> dfs = [&](vector<ui> cands, ull covered) {
    if ((covered & fullMask) == fullMask) {
      // Before recording, verify cur is not a superset of an existing solution.
      for (const auto &s : solutions)
        if (includes(cur.begin(), cur.end(), s.begin(), s.end()))
          return;
      // Remove any existing solutions that cur dominates (cur is a subset).
      solutions.erase(remove_if(solutions.begin(), solutions.end(),
                                [&](const vector<ui> &s) {
                                  return includes(s.begin(), s.end(),
                                                  cur.begin(), cur.end());
                                }),
                      solutions.end());
      solutions.push_back(cur);
      return;
    }

    const ull uncovered = fullMask & ~covered;

    // Superset pruning: cur already contains a known minimal solution so any
    // extension of cur cannot be minimal.
    for (const auto &s : solutions)
      if (includes(cur.begin(), cur.end(), s.begin(), s.end()))
        return;

    // Fail-first dead-branch check: for every uncovered constraint verify at
    // least one candidate can cover it. If any constraint is impossible, prune.
    {
      ull tmp = uncovered;
      while (tmp) {
        const int h = __builtin_ctzll(tmp);
        tmp &= tmp - 1;
        bool found = false;
        for (ui ci : cands)
          // cov is the bitmask of constraints covered by ci; check if it
          // includes h-th constraint.
          if (cov[ci] & (1ULL << h)) {
            found = true;
            break;
          }
        if (!found)
          return;
      }
    }

    for (ui ci : cands) {
      // Improvement 4: skip vertices that add no new coverage — they can
      // never appear in a minimal solution at this point.
      if (!(cov[ci] & uncovered))
        continue;

      // Build next-level candidates: those in cands with E-index > ci that
      // are adjacent to ci (enforces clique property and avoids duplicates).
      vector<ui> next;
      next.reserve(cands.size());
      for (ui cj : cands)
        if (cj > ci && compat[ci * eSize + cj])
          next.push_back(cj);

      cur.push_back(ci);
      dfs(std::move(next), covered | cov[ci]);
      cur.pop_back();
    }
  };

  dfs(std::move(initCands), forcedCov);

  // Convert E-index solutions back to actual vertex IDs, merging any forced
  // vertices that were pre-selected by unit propagation (sp1).
  vector<vector<ui>> result;
  result.reserve(solutions.size());
  for (const auto &sol : solutions) {
    vector<ui> vsol;
    vsol.reserve(sol.size() + forcedIdxs.size());
    for (ui idx : sol)
      vsol.push_back(E[idx]);
    for (ui idx : forcedIdxs)
      vsol.push_back(E[idx]);
    sort(vsol.begin(), vsol.end());
    result.push_back(std::move(vsol));
  }
  return result;
}

static string encodeClique(const vector<ui> &C) {
  return string(reinterpret_cast<const char *>(C.data()), C.size() * sizeof(ui));
}

bool ReorderSib::branchSpaceInsideClique(const vector<ui> &M,
                                         const vector<ui> &E,
                                         const vector<ui> &C) {
  return includes(C.begin(), C.end(), M.begin(), M.end()) &&
         includes(C.begin(), C.end(), E.begin(), E.end());
}

// Apply the sibling effect at this level, then continue the reorder search on
// the resulting branch seeds (mustin, expandTo).
void ReorderSib::rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo,
                       ui level, vector<char> fullSkipCheck) {
  ScopedTimer _t(rsp.rCall_ms, rsp.rCall_n);
  if (fullSkipCheck.size() != mustin.size())
    fullSkipCheck.assign(mustin.size(), 0);

  if (debug) {
    for (ui i = 0; i < level; i++)
      cout << "     ";
    cout << "Level " << level << ": mustin and expandTo sets:" << endl;
    for (ui i = 0; i < mustin.size(); i++) {
      for (ui j = 0; j < level; j++)
        cout << "     ";
      cout << "Vertex " << i << ": mustin={ ";
      for (ui v : mustin[i])
        cout << v << " ";
      cout << "}  expandTo={ ";
      for (ui v : expandTo[i])
        cout << v << " ";
      cout << "}" << endl;
    }
  }

  if (level != 0 && !expandTo.empty() && !expandTo[0].empty()) {
    vector<ui> baseM = mustin[0];
    vector<ui> baseE = expandTo[0];

    // Old cliques containing M
    vector<ui> coveringCliques = collectCoveringCliques(baseM, baseE, level);
    bool hasCoveringCliques = !coveringCliques.empty();
    bool needsFullSkip = fullSkipCheck.empty() ? false : fullSkipCheck[0];

    vector<vector<ui>> siblingSets =
        generateSiblingSetsFromCliques(baseE, coveringCliques);

    mustin.clear();
    expandTo.clear();
    fullSkipCheck.clear();
    for (const vector<ui> &S : siblingSets) {
      vector<ui> baseMustin = unionSet(baseM, S);
      vector<ui> baseExpand = commonExpand(baseE, S);
      bool branchNeedsFullSkip = needsFullSkip || hasCoveringCliques;

      if (!hasCoveringCliques || baseExpand.empty()) {
        mustin.push_back(baseMustin);
        expandTo.push_back(baseExpand);
        fullSkipCheck.push_back(branchNeedsFullSkip);
      } else {
        for (ui v : baseExpand) {
          // Direct sorted insert of v into baseMustin — avoids temporary {v}.
          vector<ui> childMustin = baseMustin;
          childMustin.insert(
              lower_bound(childMustin.begin(), childMustin.end(), v), v);
          mustin.push_back(std::move(childMustin));
          expandTo.push_back(intersect(baseExpand, adjList[v]));
          fullSkipCheck.push_back(branchNeedsFullSkip);
        }
      }
    }

    // Deduplicate branches: two sibling paths can reach the same mustin from
    // different directions (e.g. hitSet {5,9} expands as "5 then 9" and
    // "9 then 5").  Only possible when covering cliques exist; singleton
    // branches always have distinct mustins so the dedup is skipped.
    if (hasCoveringCliques) {
      unordered_map<string, ui> mustinIndex;
      vector<vector<ui>> dedupMustin;
      vector<vector<ui>> dedupExpand;
      vector<char> dedupSkip;
      for (ui i = 0; i < (ui)mustin.size(); i++) {
        string key = encodeClique(mustin[i]);
        auto [it, inserted] = mustinIndex.emplace(key, (ui)dedupMustin.size());
        if (inserted) {
          dedupMustin.push_back(std::move(mustin[i]));
          dedupExpand.push_back(std::move(expandTo[i]));
          dedupSkip.push_back(fullSkipCheck[i]);
        } else {
          // merge expandTo so no candidate is dropped
          dedupExpand[it->second] =
              unionSet(dedupExpand[it->second], expandTo[i]);
        }
      }
      mustin = std::move(dedupMustin);
      expandTo = std::move(dedupExpand);
      fullSkipCheck = std::move(dedupSkip);
    }

    if (debug) {
      for (ui i = 0; i < level; i++)
        cout << "     ";
      cout << "Level " << level << ": After Sibling Effect:" << endl;
      for (ui i = 0; i < mustin.size(); i++) {
        for (ui j = 0; j < level; j++)
          cout << "     ";
        cout << "Vertex " << i << ": mustin={ ";
        for (ui v : mustin[i])
          cout << v << " ";
        cout << "}  expandTo={ ";
        for (ui v : expandTo[i])
          cout << v << " ";
        cout << "}" << endl;
      }
    }
  }

  // Enumerate each branch until we find the next maximal clique.
  for (ui i = 0; i < (ui)mustin.size(); i++) {
    // sp5: skip branches that can't possibly produce a clique of size > 2.
    if (sp5 && mustin[i].size() + expandTo[i].size() <= 2)
      continue;
    // Empty-expandTo branches always produce exactly mustin[i] as the clique.
    // If another rCall already claimed this mustin, the result is identical —
    // skip to avoid a cross-invocation duplicate that branch-dedup can't catch.
    if (expandTo[i].empty()) {
      string key = encodeClique(mustin[i]);
      if (!claimedEmptyBranches.insert(key).second) {
        dupBlocked++;
        continue;
      }
    }
    vector<ui> R = mustin[i];
    vector<ui> Q = expandTo[i];
    bool done = false;
    enumerate(R, Q, {}, mustin, expandTo, fullSkipCheck, i, level, done);
    if (done)
      break;
  }
}

// adding candidates to Partial solution until no candidates are left.
// once maximal clique is found, reorder the remaining branches.
// X: already-processed vertices adjacent to all of R (exclusion set for maximality).
void ReorderSib::enumerate(vector<ui> &R, const vector<ui> &P,
                           const vector<ui> &X,
                           vector<vector<ui>> &mustin,
                           vector<vector<ui>> &expandTo,
                           vector<char> &fullSkipCheck, ui treeIndex, ui level,
                           bool &done) {
  ScopedTimer _t(rsp.enumerate_ms, rsp.enumerate_n);
  const ui depth = enumDepth++;
  struct DepthGuard {
    ui &d;
    ~DepthGuard() { d--; }
  } depthGuard{enumDepth};
  checksCount++;

  if (debug) {
    for (ui i = 0; i < level; i++)
      cout << "     ";
    cout << "Level " << level << ": Checking R={ ";
    for (ui v : R)
      cout << v << " ";
    cout << "}  P={ ";
    for (ui v : P)
      cout << v << " ";
    cout << "}" << endl;
  }

  // Base case: no candidates left.
  // Report R as maximal only when X is also empty (nothing extends R).
  if (P.empty()) {
    if (X.empty() && (ui)R.size() > 2) {
      vector<ui> C = R;
      sort(C.begin(), C.end());
      claimedEmptyBranches.insert(encodeClique(C));
      cliqueCount++;
      if (debug) {
        for (ui i = 0; i < level; i++)
          cout << "   ";
        cout << "Maximal Clique Found: { ";
        for (ui v : C)
          cout << v << " ";
        cout << "}" << endl;
      }
      maxCliqueSize = max(maxCliqueSize, C.size());
      ui cliqueIdx = (ui)allCliques.size();
      allCliques.push_back(C);
      for (ui v : C) {
        if ((ui)cliquesByVertexByLevel[v].size() <= level)
          cliquesByVertexByLevel[v].resize(level + 1);
        cliquesByVertexByLevel[v][level].push_back(cliqueIdx);
        cliqueCountByVertex[v]++;
      }

      done = true;

      // Reorder Logic: use the new clique to reorder the remaining branches.
      vector<vector<ui>> newMustin;
      vector<vector<ui>> newExpandTo;
      vector<char> newFullSkipCheck;

      for (ui i = treeIndex; i < (ui)mustin.size(); i++) {
        bool skipBranch = false;
        if (i < fullSkipCheck.size() && fullSkipCheck[i]) {
          skipBranch = branchSpaceInsideClique(mustin[i], expandTo[i], C);
        } else {
          skipBranch = find(C.begin(), C.end(), mustin[i].back()) != C.end();
        }
        if (skipBranch)
          continue;

        bool usesFullSkip = i < fullSkipCheck.size() && fullSkipCheck[i];
        vector<ui> reorderedExpand;
        if (usesFullSkip) {
          reorderedExpand = setDiff(unionSet(C, expandTo[i]), mustin[i]);
          for (ui mv : mustin[i])
            reorderedExpand = intersect(reorderedExpand, adjList[mv]);
        } else {
          reorderedExpand =
              intersect(setDiff(adjList[mustin[i].back()], mustin[i]),
                        unionSet(C, expandTo[i]));
        }

        if (sp6 && reorderedExpand.empty() && mustin[i].size() <= 2)
          continue;

        newMustin.push_back(mustin[i]);
        newExpandTo.push_back(reorderedExpand);
        newFullSkipCheck.push_back(usesFullSkip);

        if (debug) {
          for (ui j = 0; j < level; j++)
            cout << "     ";
          cout << "Level " << level << ": After Reorder:" << endl;
          for (ui k = 0; k < newMustin.size(); k++) {
            for (ui j = 0; j < level; j++)
              cout << "     ";
            cout << "Vertex " << k << ": mustin={ ";
            for (ui v : newMustin[k])
              cout << v << " ";
            cout << "}  expandTo={ ";
            for (ui v : newExpandTo[k])
              cout << v << " ";
            cout << "}" << endl;
          }
        }
        rCall(std::move(newMustin), std::move(newExpandTo), level + 1,
              std::move(newFullSkipCheck));
      }
    }
    return;
  }

  // Tomita pivot from P ∪ X.
  // Adaptive scoring: walk the shorter of adjList[u] (check lab) or P (check
  // adjSet) — always O(min(deg(u), |P|)) per candidate instead of O(deg(u)).
  // Mark and restore immediately so recursive children see a clean lab state.
  for (ui v : P) lab[v] = 1;
  for (ui v : X) lab[v] = 2;
  ui pivot = P[0];
  int bestNb = -1;
  const ui pSize = (ui)P.size();
  auto pivotScore = [&](ui u) -> int {
    int nb = 0;
    if ((ui)adjList[u].size() <= pSize) {
      for (ui w : adjList[u]) if (lab[w] == 1) nb++;  // walk adj, check lab
    } else {
      for (ui w : P) if (adjSet[u].count(w)) nb++;    // walk P, check adjSet
    }
    return nb;
  };
  for (ui u : P) { int nb = pivotScore(u); if (nb > bestNb) { bestNb = nb; pivot = u; } }
  for (ui u : X) { int nb = pivotScore(u); if (nb > bestNb) { bestNb = nb; pivot = u; } }
  // Restore lab before the main loop so recursive children see a clean state.
  for (ui v : P) lab[v] = 0;
  for (ui v : X) lab[v] = 0;

  // P_at_level and localX track the evolving candidate/exclusion sets.
  // intersect is used for P_new/X_new — this avoids any ancestor lab
  // interference that would arise if we used lab inside the recursive loop.
  vector<ui> &P_at_level = depthPal[depth];
  vector<ui> &localX = depthLX[depth];
  vector<ui> &P_new = depthPnew[depth];
  vector<ui> &X_new = depthXnew[depth];
  P_at_level = P;
  localX = X;
  if (localX.capacity() < X.size() + 1)
    localX.reserve(X.size() + 1);

  // Iterate exactly the branching vertices P \ N(pivot) via one merge over the
  // two sorted lists, instead of an O(1)-hash lookup for every v in P.
  const vector<ui> &pivotNbrs = adjList[pivot];
  ui nbrIdx = 0;
  for (ui v : P) {
    while (nbrIdx < pivotNbrs.size() && pivotNbrs[nbrIdx] < v)
      nbrIdx++;
    if (nbrIdx < pivotNbrs.size() && pivotNbrs[nbrIdx] == v)
      continue;

    intersectInto(P_new, P_at_level, adjList[v]);
    intersectInto(X_new, localX, adjList[v]);
    R.push_back(v);
    enumerate(R, P_new, X_new, mustin, expandTo, fullSkipCheck, treeIndex,
              level, done);
    R.pop_back();
    P_at_level.erase(lower_bound(P_at_level.begin(), P_at_level.end(), v));
    localX.insert(lower_bound(localX.begin(), localX.end(), v), v);
    if (done)
      return;
  }
}

void ReorderSib::findAllMaximalCliques() {
  rsp.reset();
  cliqueCount = 0;
  dupBlocked = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  allCliques.clear();
  claimedEmptyBranches.clear();
  cliquesByVertexByLevel.assign(n, {});
  fill(cliqueCountByVertex.begin(), cliqueCountByVertex.end(), 0);
  enumDepth = 0;

  // Tree Must-In and Expand-To sets for the initial call
  // Mustin : already a part of the clique.
  vector<vector<ui>> mustin;
  // ExpandTo : candidates that can be added to the clique.
  vector<vector<ui>> expandTo;
  for (ui v = 0; v < n; v++) {
    mustin.push_back({v});
    expandTo.push_back(adjList2[v]);
  }

  auto t0 = chrono::high_resolution_clock::now();
  // how to check if current branch is fully inside a clique:
  // fullSkipCheck[i] = flase we check using last vertex of mustin[i] in C,
  // (normal branches),  true we check using whole mustin[i] and expandTo[i] in
  // C. (sibling brnahces)
  vector<char> fullSkipCheck(n, 0);
  rCall(std::move(mustin), std::move(expandTo), 0, std::move(fullSkipCheck));
  auto t1 = chrono::high_resolution_clock::now();
  double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << fixed << setprecision(3) << "ReorderSib: cliques=" << cliqueCount
       << "  dups=" << dupBlocked << "  maxSize=" << maxCliqueSize
       << "  checks=" << checksCount << "  time=" << ms << " ms" << endl;
#if PROFILING
  rsp.print(ms);
#endif
}
