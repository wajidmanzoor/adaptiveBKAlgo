#include "../inc/helpers.h"
#include <chrono>
#include <functional>
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
// Returns true if P is a subset of any already-found maximal clique
bool PivotBK::isPSubsetOfFoundClique(const vector<ui> &P) {
  for (const vector<ui> &clique : foundCliques) {
    // Check if every vertex in P exists in this clique
    bool allFound = true;
    for (ui v : P) {
      if (find(clique.begin(), clique.end(), v) == clique.end()) {
        allFound = false;
        break;
      }
    }
    if (allFound)
      return true;
  }
  return false;
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

  // Track redundent checks for profiling
#if debug
  if (isPSubsetOfFoundClique(P)) {
    redendantChecks.push_back(P);
    redundancy++;
  }
#endif

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

  // keeps track of redunant checks (check that are subset of cliques arleady
  // found) for profiling
  redundancy = 0;
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
  cout << "Time: " << ms << " ms" << endl;

#if debug
  cout << endl << "Redundancy Count: " << redundancy << endl;

  cout << "Redundant Checks: ";
  for (const auto &check : redendantChecks) {
    cout << "{ ";
    for (ui v : check)
      cout << v << " ";
    cout << "} ";
  }
  cout << endl;
#endif
}

// ── ReorderSib profiling
// ──────────────────────────────────────────────────────
struct RSibProf {
  double rCall_ms = 0, enumerate_ms = 0, collect_ms = 0, solver_ms = 0;
  double minimal_ms = 0, commonExp_ms = 0, dedup_ms = 0, buildHit_ms = 0;
  long rCall_n = 0, enumerate_n = 0, collect_n = 0, solver_n = 0;
  long minimal_n = 0, commonExp_n = 0, dedup_n = 0, buildHit_n = 0;
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
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "solver",
           solver_ms, pct(solver_ms), solver_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "collectCoveringCliques",
           collect_ms, pct(collect_ms), collect_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "buildHitSets",
           buildHit_ms, pct(buildHit_ms), buildHit_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "minimalByInclusion",
           minimal_ms, pct(minimal_ms), minimal_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "commonExpand",
           commonExp_ms, pct(commonExp_ms), commonExp_n);
    printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", "seenCliques dedup",
           dedup_ms, pct(dedup_ms), dedup_n);
    printf("  %-30s %9.3f ms  %5.1f%%\n", "TOTAL (wall)", total_ms, 100.0);
    printf("──────────────────────────────────────────────────────────────\n");
  }
};
static RSibProf rsp;

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

// ReorderSib Implementation
ReorderSib::ReorderSib(Graph &g, DegOrder order, SibMethod method,
                       ui hitSetLimit)
    : hitSetLimit(hitSetLimit) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  // methods to get the sibling effect.
  this->method = method;
  cliquesByVertex.resize(n);

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

vector<ui> ReorderSib::compliment(const vector<ui> &vector1) {
  vector<char> seen(n, 0);
  for (ui x : vector1)
    seen[x] = 1;
  vector<ui> C;
  for (ui i = 0; i < n; i++)
    if (!seen[i])
      C.push_back(i);
  return C;
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
  for (ui v : S)
    result = intersect(result, adjList[v]);
  return result;
}

// Find previously discovered cliques that already contain M. These are the
// only cliques that matter for the sibling effect at this branch.
vector<ui> ReorderSib::collectCoveringCliques(const vector<ui> &M, ui level) {
  ScopedTimer _t(rsp.collect_ms, rsp.collect_n);
  vector<ui> result;
  if (M.empty())
    return result;

  ui seed = M[0];
  for (ui v : M) {
    if (cliquesByVertex[v].size() < cliquesByVertex[seed].size())
      seed = v;
  }

  for (ui cId : cliquesByVertex[seed]) {
    // Only keep cliques from the previous level or deeper; older levels have
    // already contributed their reorder/sibling information.
    if (cId >= allCliques.size() || foundLevel[cId] < level - 1)
      continue;

    bool containsAllMustin = includes(
        allCliques[cId].begin(), allCliques[cId].end(), M.begin(), M.end());
    if (containsAllMustin)
      result.push_back(cId);
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
  // Tighter constraints → solver has fewer candidates per constraint → runs faster.
  const vector<ui> *ids = &cliqueIds;
  vector<ui> sorted;
  if (cliqueIds.size() > maxHitSets) {
    sorted = cliqueIds;
    // Sort by clique size descending: larger clique → smaller E\C → tighter constraint
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

// Convert the sibling effect into a clique-constrained hitting-set problem
// over the already discovered covering cliques, then solve it with the
// requested strategy.
vector<vector<ui>>
ReorderSib::generateSiblingSetsFromCliques(const vector<ui> &E,
                                           const vector<ui> &cliqueIds) {
  if (cliqueIds.empty())
    return singletonBranches(E);

  vector<vector<ui>> hitSets = buildHitSets(E, cliqueIds, hitSetLimit);

  // An empty hit set means some old clique fully contains E, so no sibling
  // choice can separate the current branch from that clique.
  for (const vector<ui> &hitSet : hitSets) {
    if (hitSet.empty())
      return {};
  }

  if (cliqueIds.size() == 1)
    return singletonBranches(hitSets[0]);

  switch (method) {
  case SibMethod::BRUTE_FORCE:
    return bruteForceBySize(E, hitSets);
  case SibMethod::BACKTRACKING:
    return backtrackingBranchBound(E, hitSets);
  case SibMethod::GREEDY:
    return greedyApproximation(E, hitSets);
  case SibMethod::BITMASK:
    return bitmaskExactSearch(E, hitSets);
  case SibMethod::MIN_HITTING_SET:
    return minimumCliqueHittingSet(E, hitSets);
  case SibMethod::OPTIMIZED:
    return efficientHittingSet(E, hitSets);
  }

  return bruteForceBySize(E, hitSets);
}

vector<vector<ui>>
ReorderSib::bruteForceBySize(const vector<ui> &E,
                             const vector<vector<ui>> &hitSets) {
  vector<vector<ui>> solutions;
  vector<ui> current;

  // Enumerate clique-compatible subsets of E in increasing size.
  function<void(ui, ui)> choose = [&](ui start, ui remaining) {
    if (remaining == 0) {
      if (hitsAll(current, hitSets))
        solutions.push_back(current);
      return;
    }
    if (E.size() - start < remaining)
      return;

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
      choose(i + 1, remaining - 1);
      current.pop_back();
    }
  };

  for (ui targetSize = 1; targetSize <= E.size(); targetSize++) {
    choose(0, targetSize);
  }
  return minimalByInclusion(solutions);
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

vector<vector<ui>>
ReorderSib::greedyApproximation(const vector<ui> &E,
                                const vector<vector<ui>> &hitSets) {
  vector<ui> S;
  vector<char> inS(n, 0);
  vector<char> covered(hitSets.size(), 0);
  ui coveredCount = 0;

  // Greedily add the clique-compatible vertex that hits the largest number of
  // still-uncovered constraints.
  while (coveredCount < hitSets.size()) {
    ui bestVertex = numeric_limits<ui>::max();
    ui bestGain = 0;

    for (ui v : E) {
      if (inS[v])
        continue;

      bool connected = true;
      for (ui u : S) {
        if (!binary_search(adjList[u].begin(), adjList[u].end(), v)) {
          connected = false;
          break;
        }
      }
      if (!connected)
        continue;

      ui gain = 0;
      for (ui i = 0; i < hitSets.size(); i++) {
        if (!covered[i] &&
            binary_search(hitSets[i].begin(), hitSets[i].end(), v))
          gain++;
      }
      if (gain > bestGain || (gain == bestGain && gain > 0 && v < bestVertex)) {
        bestGain = gain;
        bestVertex = v;
      }
    }

    if (bestGain == 0)
      return {};

    S.push_back(bestVertex);
    inS[bestVertex] = 1;
    for (ui i = 0; i < hitSets.size(); i++) {
      if (!covered[i] &&
          binary_search(hitSets[i].begin(), hitSets[i].end(), bestVertex)) {
        covered[i] = 1;
        coveredCount++;
      }
    }
  }

  sort(S.begin(), S.end());
  return {S};
}

vector<vector<ui>>
ReorderSib::bitmaskExactSearch(const vector<ui> &E,
                               const vector<vector<ui>> &hitSets) {
  if (hitSets.size() > 63)
    return backtrackingBranchBound(E, hitSets);

  // Encode "which hit sets does this vertex cover?" as a 64-bit mask so the
  // exact search can update coverage with fast bitwise OR operations.
  ull fullMask = 0;
  for (ui i = 0; i < hitSets.size(); i++)
    fullMask |= (1ULL << i);

  vector<ull> masks(E.size(), 0);
  vector<ui> eIndex(n, n);
  for (ui i = 0; i < E.size(); i++)
    eIndex[E[i]] = i;
  for (ui j = 0; j < hitSets.size(); j++) {
    for (ui v : hitSets[j]) {
      ui idx = eIndex[v];
      if (idx < E.size())
        masks[idx] |= (1ULL << j);
    }
  }

  vector<vector<ui>> solutions;
  vector<ui> current;
  vector<ui> currentIndices;
  const size_t eSize = E.size();
  const bool useCompatTable = eSize <= 4096;
  vector<char> compatible;
  if (useCompatTable) {
    compatible.assign(eSize * eSize, 0);
    for (ui i = 0; i < eSize; i++) {
      for (ui j = 0; j < i; j++) {
        if (binary_search(adjList[E[j]].begin(), adjList[E[j]].end(), E[i])) {
          compatible[i * eSize + j] = 1;
          compatible[j * eSize + i] = 1;
        }
      }
    }
  }

  function<void(ui, ui, ull)> dfs = [&](ui start, ui remaining, ull mask) {
    if (remaining == 0) {
      if (mask == fullMask)
        solutions.push_back(current);
      return;
    }
    if (E.size() - start < remaining)
      return;

    for (ui i = start; i < E.size(); i++) {
      bool connected = true;
      if (useCompatTable) {
        for (ui idx : currentIndices) {
          if (!compatible[idx * eSize + i]) {
            connected = false;
            break;
          }
        }
      } else {
        for (ui v : current) {
          if (!binary_search(adjList[v].begin(), adjList[v].end(), E[i])) {
            connected = false;
            break;
          }
        }
      }
      if (!connected)
        continue;

      current.push_back(E[i]);
      currentIndices.push_back(i);
      dfs(i + 1, remaining - 1, mask | masks[i]);
      currentIndices.pop_back();
      current.pop_back();
    }
  };

  for (ui targetSize = 1; targetSize <= E.size(); targetSize++) {
    dfs(0, targetSize, 0);
  }
  return minimalByInclusion(solutions);
}

vector<vector<ui>>
ReorderSib::minimumCliqueHittingSet(const vector<ui> &E,
                                    const vector<vector<ui>> &hitSets) {
  return backtrackingBranchBound(E, hitSets);
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
  const ull fullMask = (hSize == 64) ? ~0ULL : ((1ULL << hSize) - 1);
  vector<ull> cov(eSize, 0);
  {
    vector<ui> eIdx(n, eSize);
    for (ui i = 0; i < eSize; i++)
      eIdx[E[i]] = i;
    for (ui h = 0; h < hSize; h++)
      for (ui v : hitSets[h])
        if (eIdx[v] < eSize)
          cov[eIdx[v]] |= (1ULL << h);
  }

  // Initial candidate order: descending coverage count so high-utility
  // vertices are tried first, finding solutions sooner for better pruning.
  vector<ui> initCands(eSize);
  iota(initCands.begin(), initCands.end(), 0);
  sort(initCands.begin(), initCands.end(), [&](ui a, ui b) {
    return __builtin_popcountll(cov[a]) > __builtin_popcountll(cov[b]);
  });

  // ── DFS ───────────────────────────────────────────────────────────────────

  // solutions is maintained as a live minimal-by-inclusion set throughout.
  vector<vector<ui>> solutions;
  // cur holds E-indices of the partial solution, in strictly increasing order.
  vector<ui> cur;

  // cands : E-indices still reachable (clique-compatible with cur, index >
  //         last element of cur). Passed by value so each level owns its copy.
  // covered: bitmask of constraints already satisfied by cur.
  function<void(vector<ui>, ull)> dfs = [&](vector<ui> cands, ull covered) {
    if (covered == fullMask) {
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

  dfs(std::move(initCands), 0);

  // Convert E-index solutions back to actual vertex IDs.
  vector<vector<ui>> result;
  result.reserve(solutions.size());
  for (const auto &sol : solutions) {
    vector<ui> vsol;
    vsol.reserve(sol.size());
    for (ui idx : sol)
      vsol.push_back(E[idx]);
    sort(vsol.begin(), vsol.end());
    result.push_back(std::move(vsol));
  }
  return result;
}

static string encodeClique(const vector<ui> &C) {
  string key;
  key.reserve(C.size() * 6);
  for (ui v : C) {
    key += to_string(v);
    key.push_back(',');
  }
  return key;
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

    // Old cliques containing M are exactly the cliques that can force a
    // sibling split at this branch.
    vector<ui> coveringCliques = collectCoveringCliques(baseM, level);
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
          vector<ui> childMustin = unionSet(baseMustin, {v});
          mustin.push_back(childMustin);
          expandTo.push_back(intersect(baseExpand, adjList[v]));
          fullSkipCheck.push_back(branchNeedsFullSkip);
        }
      }
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

  // Enumerate each  branch until we find the next maximal clique.
  for (ui i = 0; i < (ui)mustin.size(); i++) {
    vector<ui> R = mustin[i];
    vector<ui> Q = expandTo[i];
    bool done = false;
    enumerate(R, Q, mustin, expandTo, fullSkipCheck, i, level, done);
    if (done)
      break;
  }
}

// adding candidates to Partial solution untill no candidates are left.
// once maximal clique is found, reorder the remaining branches.
void ReorderSib::enumerate(vector<ui> &R, vector<ui> &Q,
                           vector<vector<ui>> &mustin,
                           vector<vector<ui>> &expandTo,
                           vector<char> &fullSkipCheck, ui treeIndex, ui level,
                           bool &done) {
  ScopedTimer _t(rsp.enumerate_ms, rsp.enumerate_n);
  checksCount++;

  if (debug) {
    for (ui i = 0; i < level; i++)
      cout << "     ";
    cout << "Level " << level << ": Checking R={ ";
    for (ui v : R)
      cout << v << " ";
    cout << "}  Q={ ";
    for (ui v : Q)
      cout << v << " ";
    cout << "}" << endl;
  }

  // is no candidate left to expand, R is maximal clique.
  if (Q.empty()) {
    if ((ui)R.size() > 2) {
      vector<ui> C = R;
      sort(C.begin(), C.end());
      { ScopedTimer _td(rsp.dedup_ms, rsp.dedup_n);
        if (!seenCliques.insert(encodeClique(C)).second) return; }
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
      foundLevel.push_back(level);
      for (ui v : C)
        cliquesByVertex[v].push_back(cliqueIdx);

      // stops the DFS in enumeerate and tree expansion in rCall.
      done = true;

      // Reorder Logic: use the new clique to reorder the remaining branches in
      // this tree.
      vector<vector<ui>> newMustin;
      vector<vector<ui>> newExpandTo;
      vector<char> newFullSkipCheck;

      for (ui i = treeIndex; i < (ui)mustin.size(); i++) {
        bool skipBranch = false;
        if (i < fullSkipCheck.size() && fullSkipCheck[i]) {
          // Full-skip mode means the whole seeded search space for this branch
          // sits inside the clique we just found, so we can discard it.
          skipBranch = branchSpaceInsideClique(mustin[i], expandTo[i], C);
        } else {
          skipBranch = find(C.begin(), C.end(), mustin[i].back()) != C.end();
        }
        if (skipBranch)
          continue;

        bool usesFullSkip = i < fullSkipCheck.size() && fullSkipCheck[i];
        vector<ui> reorderedExpand;
        if (usesFullSkip) {
          // Rebuild the entire branch space from M and the new clique C.
          reorderedExpand = setDiff(unionSet(C, expandTo[i]), mustin[i]);
          for (ui mv : mustin[i])
            reorderedExpand = intersect(reorderedExpand, adjList[mv]);
        } else {
          // Base reorder update: branch from the last mandatory vertex only.
          reorderedExpand =
              intersect(setDiff(adjList[mustin[i].back()], mustin[i]),
                        unionSet(C, expandTo[i]));
        }

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
      return;
    }
  }

  // Enemerate candidates
  for (ui v : Q) {
    R.push_back(v);
    vector<ui> Qp = intersect(Q, adjList2[v]);
    enumerate(R, Qp, mustin, expandTo, fullSkipCheck, treeIndex, level, done);
    R.pop_back();
    if (done)
      return;
  }
}

void ReorderSib::findAllMaximalCliques() {
  rsp.reset();
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  allCliques.clear();
  foundLevel.clear();
  seenCliques.clear();
  for (ui v = 0; v < n; v++)
    cliquesByVertex[v].clear();

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
  vector<char> fullSkipCheck(n, 0);
  rCall(std::move(mustin), std::move(expandTo), 0, std::move(fullSkipCheck));
  auto t1 = chrono::high_resolution_clock::now();
  double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << "ReorderSib: cliques=" << cliqueCount << "  maxSize=" << maxCliqueSize
       << "  checks=" << checksCount << "  time=" << ms << " ms" << endl;
  rsp.print(ms);
}
