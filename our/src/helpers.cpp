#include "../inc/helpers.h"
#include "../inc/fast_plex3.h"
#include <chrono>
#include <functional>
#include <iomanip>
#include <numeric>
// Returns peelSeq index of verticies by core value
// peelSeq[0] = highest-core vertex, peelSeq[n-1] = lowest.
static vector<ui> computePeelSeq(const Graph &g) {
  ui n = g.n;
  if (n == 0)
    return {};

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
  if (P.empty() && X.empty()) {
    if (R.size() <= 2)
      return;

    // Found a maximal clique
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());
    return;
  }

  // P empty but X non-empty: R is NOT maximal, prune
  if (P.empty())
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
  double collect_ms = 0, solver_ms = 0;
  double commonExp_ms = 0, buildHit_ms = 0;
  double intersect_ms = 0, setdiff_ms = 0, unionset_ms = 0, encode_ms = 0;
  long collect_n = 0, solver_n = 0;
  long commonExp_n = 0, buildHit_n = 0;
  long intersect_n = 0, setdiff_n = 0, unionset_n = 0, encode_n = 0;
  ull solver_hsize_sum = 0, solver_esize_sum = 0;
  ull solver_compat_eligible = 0, solver_compat_survivors = 0;
  long solver_h_gt63 = 0, solver_h_gt128 = 0;
  long solver_h_le8 = 0, solver_h_le16 = 0, solver_h_le32 = 0;
  long solver_h_le63 = 0, solver_h_le128 = 0, solver_h_gt128_bucket = 0;
  void reset() { *this = RSibProf{}; }
  void print(double total_ms) const {
    auto pct = [&](double v) {
      return total_ms > 0 ? 100.0 * v / total_ms : 0.0;
    };
    const vector<pair<const char *, pair<double, long>>> rows = {
        {"collect covering cliques", {collect_ms, collect_n}},
        {"buildHitSets", {buildHit_ms, buildHit_n}},
        {"solver", {solver_ms, solver_n}},
        {"commonExpand", {commonExp_ms, commonExp_n}},
        {"intersect/intersectInto", {intersect_ms, intersect_n}},
        {"setDiff", {setdiff_ms, setdiff_n}},
        {"unionSet", {unionset_ms, unionset_n}},
        {"encodeClique", {encode_ms, encode_n}},
    };
    double profiled_ms = 0.0;
    for (const auto &row : rows)
      profiled_ms += row.second.first;
    const double other_ms = max(0.0, total_ms - profiled_ms);

    printf(
        "\n── ReorderSib cost breakdown ─────────────────────────────────\n");
    for (const auto &row : rows) {
      printf("  %-30s %9.3f ms  %5.1f%%  calls=%ld\n", row.first,
             row.second.first, pct(row.second.first), row.second.second);
    }
    printf("  %-30s %9.3f ms  %5.1f%%\n", "other / untimed", other_ms,
           pct(other_ms));
    printf("  %-30s %9.3f ms  %5.1f%%\n", "PROFILED subtotal", profiled_ms,
           pct(profiled_ms));
    printf("  %-30s %9.3f ms  %5.1f%%\n", "TOTAL (wall)", total_ms, 100.0);
    if (solver_n > 0) {
      const double avgE = static_cast<double>(solver_esize_sum) / solver_n;
      const double avgH = static_cast<double>(solver_hsize_sum) / solver_n;
      const double compatPct =
          solver_compat_eligible > 0
              ? 100.0 * static_cast<double>(solver_compat_survivors) /
                    solver_compat_eligible
              : 0.0;
      printf("\n  Solver stats\n");
      printf("  %-30s %.2f\n", "avg |E|", avgE);
      printf("  %-30s %.2f\n", "avg hSize", avgH);
      printf("  %-30s %ld\n", "calls with hSize > 63", solver_h_gt63);
      printf("  %-30s %ld\n", "calls with hSize > 128", solver_h_gt128);
      printf("  %-30s %llu / %llu  (%4.1f%%)\n", "compat survivors / eligible",
             solver_compat_survivors, solver_compat_eligible, compatPct);
      printf("  %-30s <=8:%ld  9-16:%ld  17-32:%ld\n", "hSize buckets",
             solver_h_le8, solver_h_le16, solver_h_le32);
      printf("  %-30s 33-63:%ld  64-128:%ld  >128:%ld\n", "", solver_h_le63,
             solver_h_le128, solver_h_gt128_bucket);
    }
    printf("──────────────────────────────────────────────────────────────\n");
  }
};
static RSibProf rsp;

#if PROFILING
struct ScopedTimer {
  using clock = chrono::high_resolution_clock;
  using duration = clock::duration;
  static thread_local ScopedTimer *current;

  chrono::high_resolution_clock::time_point t0;
  duration child = duration::zero();
  double &acc;
  long &cnt;
  ScopedTimer *parent;

  ScopedTimer(double &a, long &c)
      : t0(clock::now()), acc(a), cnt(c), parent(current) {
    current = this;
  }

  ~ScopedTimer() {
    const duration elapsed = clock::now() - t0;
    duration self = elapsed - child;
    if (self < duration::zero())
      self = duration::zero();
    acc += chrono::duration<double, milli>(self).count();
    ++cnt;
    current = parent;
    if (parent != nullptr)
      parent->child += elapsed;
  }
};
thread_local ScopedTimer *ScopedTimer::current = nullptr;
#else
struct ScopedTimer {
  ScopedTimer(double &, long &) {}
};
#endif

// ReorderSib Implementation
ReorderSib::ReorderSib(Graph &g, DegOrder order, ui minCliqueSize)
    : minCliqueSize(max<ui>(1, minCliqueSize)) {
  n = g.n;
  cliqueCount = 0;
  dupBlocked = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  solverWorkBudget = 0;
  solverBudgetFallbacks = 0;

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

  internalToOriginal.resize(n);
  for (ui original = 0; original < n; original++)
    internalToOriginal[perm[original]] = original;

  buildAdjLists(g, perm, adjList, adjList2);

  adjSet.resize(n);
  for (ui u = 0; u < n; u++)
    for (ui v : adjList[u])
      adjSet[u].insert(v);

  eIndex.assign(n, 0);
  eIndexStamp.assign(n, 0);
  eIndexToken = 0;
}

void ReorderSib::intersectInto(vector<ui> &out, const vector<ui> &A,
                               const vector<ui> &B) {
  ScopedTimer _t(rsp.intersect_ms, rsp.intersect_n);
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  if (A.size() * 8 < B.size()) {
    for (ui v : A)
      if (binary_search(B.begin(), B.end(), v))
        out.push_back(v);
    return;
  }
  if (B.size() * 8 < A.size()) {
    for (ui v : B)
      if (binary_search(A.begin(), A.end(), v))
        out.push_back(v);
    return;
  }
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
  ScopedTimer _t(rsp.setdiff_ms, rsp.setdiff_n);
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

void ReorderSib::setDiffInto(vector<ui> &out, const vector<ui> &A,
                             const vector<ui> &B) {
  ScopedTimer _t(rsp.setdiff_ms, rsp.setdiff_n);
  out.clear();
  if (out.capacity() < A.size())
    out.reserve(A.size());
  ui i = 0, j = 0;
  while (i < A.size()) {
    if (j == (ui)B.size() || A[i] < B[j]) {
      out.push_back(A[i]);
      i++;
    } else if (A[i] == B[j]) {
      i++;
      j++;
    } else {
      j++;
    }
  }
}

vector<ui> ReorderSib::unionSet(const vector<ui> &A, const vector<ui> &B) {
  ScopedTimer _t(rsp.unionset_ms, rsp.unionset_n);
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

void ReorderSib::recordSolverCallStats(ui eSize, ui hSize) {
  rsp.solver_esize_sum += eSize;
  rsp.solver_hsize_sum += hSize;
  if (hSize > 63)
    rsp.solver_h_gt63++;
  if (hSize > 128)
    rsp.solver_h_gt128++;

  if (hSize <= 8)
    rsp.solver_h_le8++;
  else if (hSize <= 16)
    rsp.solver_h_le16++;
  else if (hSize <= 32)
    rsp.solver_h_le32++;
  else if (hSize <= 63)
    rsp.solver_h_le63++;
  else if (hSize <= 128)
    rsp.solver_h_le128++;
  else
    rsp.solver_h_gt128_bucket++;
}

void ReorderSib::recordSolverCompatStats(ull eligible, ull survivors) {
  rsp.solver_compat_eligible += eligible;
  rsp.solver_compat_survivors += survivors;
}

vector<ui> ReorderSib::commonExpand(const vector<ui> &E, const vector<ui> &S) {
  ScopedTimer _t(rsp.commonExp_ms, rsp.commonExp_n);
  if (S.empty())
    return E;
  if (S.size() == 1) {
    vector<ui> result;
    intersectInto(result, E, adjList[S[0]]);
    return result;
  }

  vector<ui> order = S;
  sort(order.begin(), order.end(),
       [&](ui a, ui b) { return adjList[a].size() < adjList[b].size(); });

  vector<ui> result;
  intersectInto(result, E, adjList[order[0]]);
  vector<ui> filtered;
  setDiffInto(filtered, result, S);
  result.swap(filtered);

  vector<ui> scratch;
  for (ui idx = 1; idx < (ui)order.size() && !result.empty(); idx++) {
    intersectInto(scratch, result, adjList[order[idx]]);
    result.swap(scratch);
  }
  return result;
}

// Proof-faithful cover lookup for Pure-ReorderSib.  Unlike the legacy
// collector, this path deliberately ignores the recent-level and weak-cover
// filters: every emitted clique containing M must be visible before FindOne
// may run.
vector<ui> ReorderSib::collectAllCoveringCliques(const vector<ui> &M) {
  ScopedTimer _t(rsp.collect_ms, rsp.collect_n);
  vector<ui> result;
  if (M.empty())
    return result;

  ui seed = M[0];
  for (ui v : M)
    if (cliqueCountByVertex[v] < cliqueCountByVertex[seed])
      seed = v;

  const auto &buckets = cliquesByVertexByLevel[seed];
  if (cliqueCountByVertex[seed] >
      static_cast<ull>(numeric_limits<size_t>::max()))
    throw overflow_error("per-vertex clique index exceeds size_t");
  result.reserve(static_cast<size_t>(cliqueCountByVertex[seed]));
  for (const auto &bucket : buckets) {
    for (ui cliqueId : bucket) {
      const vector<ui> &C = allCliques[cliqueId];
      if (includes(C.begin(), C.end(), M.begin(), M.end()))
        result.push_back(cliqueId);
    }
  }
  return result;
}

// Each covering clique C contributes the constraint "pick something from E
// that is outside C", which is exactly E \ C.
vector<vector<ui>> ReorderSib::buildHitSets(const vector<ui> &E,
                                            const vector<ui> &cliqueIds) {
  ScopedTimer _t(rsp.buildHit_ms, rsp.buildHit_n);

  vector<vector<ui>> hitSets;
  hitSets.reserve(cliqueIds.size());
  for (ui cId : cliqueIds)
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

// Exact sibling-seed generation used by the theorem-aligned Pure lane.  It
// applies every covering-clique constraint.
vector<vector<ui>>
ReorderSib::generateExactSiblingSets(const vector<ui> &E,
                                     const vector<ui> &coveringCliqueIds,
                                     bool *usePivotFallback) {
  if (usePivotFallback != nullptr)
    *usePivotFallback = false;
  if (coveringCliqueIds.empty())
    return {};

  vector<vector<ui>> hitSets = buildHitSets(E, coveringCliqueIds);
  for (const vector<ui> &hitSet : hitSets)
    if (hitSet.empty())
      return {};

  if (coveringCliqueIds.size() == 1)
    return singletonBranches(hitSets[0]);
  return efficientHittingSet(E, hitSets, usePivotFallback);
}

// Optimized exact solver for all minimal clique-constrained hitting sets.
//
// Main implementation details:
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
ReorderSib::efficientHittingSet(const vector<ui> &inputE,
                                const vector<vector<ui>> &inputHitSets,
                                bool *usePivotFallback) {
  ScopedTimer _t(rsp.solver_ms, rsp.solver_n);
  if (usePivotFallback != nullptr)
    *usePivotFallback = false;

  // Preprocess on vertex sets before selecting a fixed-mask fallback.  This is
  // deliberately representation-independent so all three experiment variants
  // make their capacity decision from the same exact reduced problem.
  if (++eIndexToken == 0) {
    fill(eIndexStamp.begin(), eIndexStamp.end(), 0);
    eIndexToken = 1;
  }
  for (size_t i = 0; i < inputE.size(); ++i) {
    eIndex[inputE[i]] = static_cast<ui>(i);
    eIndexStamp[inputE[i]] = eIndexToken;
  }
  vector<char> active(inputE.size(), 1);
  vector<char> forcedMark(inputE.size(), 0);
  auto localIndexOf = [&](ui v) -> size_t {
    return static_cast<size_t>(eIndex[v]);
  };

  vector<vector<ui>> hitSets = inputHitSets;
  vector<ui> preForced;
  bool infeasible = false;

  auto normalizeAndSubsume = [&]() {
    for (vector<ui> &H : hitSets) {
      sort(H.begin(), H.end());
      H.erase(unique(H.begin(), H.end()), H.end());
    }
    sort(hitSets.begin(), hitSets.end());
    hitSets.erase(unique(hitSets.begin(), hitSets.end()), hitSets.end());
    sort(hitSets.begin(), hitSets.end(),
         [](const vector<ui> &a, const vector<ui> &b) {
           if (a.size() != b.size())
             return a.size() < b.size();
           return a < b;
         });

    vector<vector<ui>> kept;
    kept.reserve(hitSets.size());
    for (vector<ui> &H : hitSets) {
      bool redundant = false;
      for (const vector<ui> &tight : kept) {
        if (tight.size() > H.size())
          break;
        if (includes(H.begin(), H.end(), tight.begin(), tight.end())) {
          redundant = true;
          break;
        }
      }
      if (!redundant)
        kept.push_back(std::move(H));
    }
    hitSets = std::move(kept);
  };

  while (!infeasible) {
    vector<vector<ui>> filtered;
    filtered.reserve(hitSets.size());
    for (const vector<ui> &H : hitSets) {
      bool alreadyHit = false;
      vector<ui> remaining;
      remaining.reserve(H.size());
      for (ui v : H) {
        const size_t local = localIndexOf(v);
        if (forcedMark[local]) {
          alreadyHit = true;
          break;
        }
        if (active[local])
          remaining.push_back(v);
      }
      if (alreadyHit)
        continue;
      if (remaining.empty()) {
        infeasible = true;
        break;
      }
      filtered.push_back(std::move(remaining));
    }
    if (infeasible)
      break;

    hitSets = std::move(filtered);
    normalizeAndSubsume();

    ui forced = n;
    for (const vector<ui> &H : hitSets) {
      if (H.size() == 1) {
        forced = H[0];
        break;
      }
    }
    if (forced == n)
      break;
    const size_t forcedLocal = localIndexOf(forced);
    if (!active[forcedLocal]) {
      infeasible = true;
      break;
    }
    for (ui v : preForced) {
      if (!adj(v, forced)) {
        infeasible = true;
        break;
      }
    }
    if (infeasible)
      break;

    preForced.push_back(forced);
    forcedMark[forcedLocal] = 1;
    active[forcedLocal] = 0;
    for (size_t i = 0; i < inputE.size(); ++i)
      if (active[i] && !adj(forced, inputE[i]))
        active[i] = 0;
  }

  if (infeasible)
    return {};
  sort(preForced.begin(), preForced.end());

  if (hitSets.empty()) {
    if (preForced.empty())
      return {};
    return {preForced};
  }

  vector<char> useful(inputE.size(), 0);
  for (const vector<ui> &H : hitSets)
    for (ui v : H)
      useful[localIndexOf(v)] = 1;
  vector<ui> E;
  E.reserve(inputE.size());
  for (size_t i = 0; i < inputE.size(); ++i)
    if (active[i] && useful[i])
      E.push_back(inputE[i]);

  const ui eSize = (ui)E.size();
  const ui hSize = (ui)hitSets.size();
  recordSolverCallStats(eSize, hSize);

  static_assert(PURE_HITSET_CAPACITY == 128,
                "this minimal source keeps only the 128-bit variant");
  if (hSize > PURE_HITSET_CAPACITY) {
    if (usePivotFallback == nullptr)
      throw logic_error("128-bit overflow requires a PXR fallback target");
    *usePivotFallback = true;
    return {};
  }
  using Mask = array<ull, 2>;
  static constexpr ui maskWords = 2;
  auto zeroMask = []() -> Mask { return Mask{}; };

  auto setBit = [](Mask &m, ui bit) { m[bit >> 6] |= (1ULL << (bit & 63)); };
  auto clearBit = [](Mask &m, ui bit) { m[bit >> 6] &= ~(1ULL << (bit & 63)); };
  auto hasBit = [](const Mask &m, ui bit) -> bool {
    return (m[bit >> 6] & (1ULL << (bit & 63))) != 0;
  };
  auto anyMask = [&](const Mask &m) -> bool {
    for (ui word = 0; word < maskWords; word++)
      if (m[word] != 0)
        return true;
    return false;
  };
  auto orEq = [&](Mask &dst, const Mask &src) {
    for (ui word = 0; word < maskWords; word++)
      dst[word] |= src[word];
  };
  auto andNot = [&](const Mask &a, const Mask &b) -> Mask {
    Mask result = zeroMask();
    for (ui word = 0; word < maskWords; word++)
      result[word] = a[word] & ~b[word];
    return result;
  };
  auto intersects = [&](const Mask &a, const Mask &b) -> bool {
    for (ui word = 0; word < maskWords; word++)
      if ((a[word] & b[word]) != 0)
        return true;
    return false;
  };
  auto coversAll = [&](const Mask &covered, const Mask &need) -> bool {
    for (ui word = 0; word < maskWords; word++)
      if ((covered[word] & need[word]) != need[word])
        return false;
    return true;
  };
  auto popcountMask = [&](const Mask &m) -> int {
    int count = 0;
    for (ui word = 0; word < maskWords; word++)
      count += __builtin_popcountll(m[word]);
    return count;
  };
  auto popNextBit = [&](Mask &m) -> int {
    for (ui word = 0; word < maskWords; word++) {
      if (m[word]) {
        const int bit = __builtin_ctzll(m[word]);
        m[word] &= (m[word] - 1);
        return static_cast<int>(word * 64 + bit);
      }
    }
    return -1;
  };

  // ── Pre-computation ───────────────────────────────────────────────────────

  // compat[i*eSize+j] = 1  iff  E[i] and E[j] are adjacent in the graph.
  vector<char> compat(eSize * eSize, 0);
  for (ui i = 0; i < eSize; i++)
    for (ui j = i + 1; j < eSize; j++)
      if (binary_search(adjList[E[i]].begin(), adjList[E[i]].end(), E[j]))
        compat[i * eSize + j] = compat[j * eSize + i] = 1;

  // cov[i] = bitmask of hitSets that E[i] covers.
  // Non-const so constraint subsumption can clear implied bits.
  Mask fullMask = zeroMask();
  for (ui bit = 0; bit < hSize; bit++)
    setBit(fullMask, bit);
  vector<Mask> cov(eSize, zeroMask());
  {
    // Sort hitSets ascending by size before encoding into bitmasks.
    // Smaller hitSets are harder to satisfy (fewer candidates), so assigning
    // them to lower bit positions makes the fail-first dead-branch check
    // encounter the hardest constraint first and prune sooner.
    vector<ui> hOrder(hSize);
    iota(hOrder.begin(), hOrder.end(), 0);
    sort(hOrder.begin(), hOrder.end(),
         [&](ui a, ui b) { return hitSets[a].size() < hitSets[b].size(); });

    if (++eIndexToken == 0) {
      fill(eIndexStamp.begin(), eIndexStamp.end(), 0);
      eIndexToken = 1;
    }
    for (ui i = 0; i < eSize; i++) {
      eIndex[E[i]] = i;
      eIndexStamp[E[i]] = eIndexToken;
    }
    for (ui bit = 0; bit < hSize; bit++) {
      ui h = hOrder[bit];
      for (ui v : hitSets[h]) {
        if (eIndexStamp[v] == eIndexToken)
          setBit(cov[eIndex[v]], bit);
      }
    }
  }

  // Initial candidate order: descending coverage count so high-utility
  // vertices are tried first, finding solutions sooner for better pruning.
  vector<ui> initCands(eSize);
  iota(initCands.begin(), initCands.end(), 0);
  sort(initCands.begin(), initCands.end(),
       [&](ui a, ui b) { return popcountMask(cov[a]) > popcountMask(cov[b]); });

  // Unit propagation.
  // For any constraint covered by exactly one candidate, that candidate is
  // forced into every solution. Pre-select all forced candidates, filter the
  // remaining candidates to be clique-compatible with them, and start the DFS
  // with the forced coverage already accumulated.
  vector<ui> forcedIdxs; // E-indices forced into every solution
  Mask forcedCov = zeroMask();

  {
    vector<ui> activeCands = initCands;
    bool changed = true;
    bool conflict = false;

    while (changed && !conflict) {
      changed = false;
      Mask uncov = andNot(fullMask, forcedCov);
      if (!anyMask(uncov))
        break;

      Mask tmp = uncov;
      while (!conflict) {
        const int h = popNextBit(tmp);
        if (h < 0)
          break;

        ui sole = eSize;
        int cnt = 0;
        for (ui ci : activeCands) {
          if (hasBit(cov[ci], (ui)h)) {
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
          orEq(forcedCov, cov[sole]);

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

  // Constraint subsumption.
  // Constraint i is subsumed by constraint j when every candidate covering j
  // also covers i (coverSet(j) ⊆ coverSet(i)). Satisfying j then implies
  // satisfying i, so i can be dropped from fullMask.
  // Work on the post-propagation candidates so forced coverage is reflected.
  {
    for (ui h = 0; h < hSize; h++) {
      if (!hasBit(fullMask, h))
        continue; // already dropped
      if (hasBit(forcedCov, h))
        continue; // h already covered by forced, DFS won't see it
      for (ui g = 0; g < hSize; g++) {
        if (g == h || !hasBit(fullMask, g))
          continue;
        // g must not be force-covered: if no initCand covers g (because a
        // forced vertex was the sole cover and propagation removed it), the
        // check would pass vacuously and incorrectly drop h.
        if (hasBit(forcedCov, g))
          continue;
        // Check coverSet(g, initCands) ⊆ coverSet(h, initCands):
        // no remaining candidate covers g but not h.
        bool gSubsumesH = true;
        for (ui ci : initCands) {
          if (hasBit(cov[ci], g) && !hasBit(cov[ci], h)) {
            gSubsumesH = false;
            break;
          }
        }
        if (gSubsumesH) {
          clearBit(fullMask, h); // drop h — implied by g
          break;
        }
      }
    }
  }

  // Early exit when forced coverage satisfies all remaining constraints.
  if (coversAll(forcedCov, fullMask)) {
    if (forcedIdxs.empty() && preForced.empty())
      return {}; // nothing forced, nothing to cover
    vector<ui> sol = preForced;
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
  vector<vector<ui>> candScratch(eSize + 1);
  ull solverWork = 0;
  bool budgetExceeded = false;

  // cands : E-indices still reachable (clique-compatible with cur, index >
  //         last element of cur).
  // covered: bitmask of constraints already satisfied by cur.
  function<void(const vector<ui> &, Mask, ui)> dfs =
      [&](const vector<ui> &cands, Mask covered, ui depth) {
        if (budgetExceeded)
          return;
        if (coversAll(covered, fullMask)) {
          // Before recording, verify cur is not a superset of an existing
          // solution.
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

        const Mask uncovered = andNot(fullMask, covered);

        // Superset pruning: cur already contains a known minimal solution so
        // any extension of cur cannot be minimal.
        for (const auto &s : solutions)
          if (includes(cur.begin(), cur.end(), s.begin(), s.end()))
            return;

        // Fail-first dead-branch check: for every uncovered constraint verify
        // at least one candidate can cover it. If any constraint is impossible,
        // prune.
        {
          Mask tmp = uncovered;
          while (true) {
            const int h = popNextBit(tmp);
            if (h < 0)
              break;
            bool found = false;
            for (ui ci : cands)
              // cov is the bitmask of constraints covered by ci; check if it
              // includes h-th constraint.
              if (hasBit(cov[ci], (ui)h)) {
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
          if (!intersects(cov[ci], uncovered))
            continue;

          // Build next-level candidates: those in cands with E-index > ci that
          // are adjacent to ci (enforces clique property and avoids
          // duplicates).
          vector<ui> &next = candScratch[depth];
          next.clear();
          if (next.capacity() < cands.size())
            next.reserve(cands.size());
          ull eligible = 0;
          for (ui cj : cands)
            if (cj > ci) {
              if (solverWorkBudget != 0 && usePivotFallback != nullptr &&
                  solverWork >= solverWorkBudget) {
                budgetExceeded = true;
                break;
              }
              ++solverWork;
              eligible++;
              if (compat[ci * eSize + cj])
                next.push_back(cj);
            }
          if (budgetExceeded)
            return;
          recordSolverCompatStats(eligible, (ull)next.size());

          cur.push_back(ci);
          Mask nextCovered = covered;
          orEq(nextCovered, cov[ci]);
          dfs(next, nextCovered, depth + 1);
          cur.pop_back();
          if (budgetExceeded)
            return;
        }
      };

  dfs(initCands, forcedCov, 1);

  if (budgetExceeded) {
    *usePivotFallback = true;
    addCliqueCountOrThrow(solverBudgetFallbacks, 1);
    return {};
  }

  // Convert E-index solutions back to actual vertex IDs, merging any forced
  // vertices that were pre-selected by unit propagation.
  vector<vector<ui>> result;
  result.reserve(solutions.size());
  for (const auto &sol : solutions) {
    vector<ui> vsol;
    vsol.reserve(sol.size() + forcedIdxs.size() + preForced.size());
    for (ui idx : sol)
      vsol.push_back(E[idx]);
    for (ui idx : forcedIdxs)
      vsol.push_back(E[idx]);
    vsol.insert(vsol.end(), preForced.begin(), preForced.end());
    sort(vsol.begin(), vsol.end());
    result.push_back(std::move(vsol));
  }
  return result;
}

ui ReorderSib::pureNeighborsInP(ui u, const vector<ui> &P) const {
  ui score = 0;
  if (adjList[u].size() < P.size()) {
    for (ui v : adjList[u])
      score += binary_search(P.begin(), P.end(), v);
  } else {
    for (ui v : P)
      score += adjSet[u].count(v) != 0;
  }
  return score;
}

void ReorderSib::scanPurePXRState(const vector<ui> &P, const vector<ui> &X,
                                  ui &pivot, ui &minPScore, ui &universalP,
                                  bool &xUniversal) const {
  const ui pSize = static_cast<ui>(P.size());
  pivot = P.front();
  minPScore = pSize;
  universalP = numeric_limits<ui>::max();
  xUniversal = false;
  int bestScore = -1;

  for (ui u : P) {
    const ui score = pureNeighborsInP(u, P);
    minPScore = min(minPScore, score);
    if (score + 1 == pSize && universalP == numeric_limits<ui>::max())
      universalP = u;
    if (static_cast<int>(score) > bestScore) {
      bestScore = static_cast<int>(score);
      pivot = u;
    }
  }
  for (ui u : X) {
    const ui score = pureNeighborsInP(u, P);
    xUniversal = xUniversal || score == pSize;
    if (static_cast<int>(score) > bestScore) {
      bestScore = static_cast<int>(score);
      pivot = u;
    }
  }
}

void ReorderSib::pureMatchingParts(const vector<ui> &P, vector<ui> &forced,
                                   vector<pair<ui, ui>> &missingEdges) const {
  forced.clear();
  missingEdges.clear();
  vector<char> paired(P.size(), 0);
  for (size_t i = 0; i < P.size(); ++i) {
    if (paired[i])
      continue;
    size_t mate = P.size();
    for (size_t j = i + 1; j < P.size(); ++j) {
      if (!paired[j] && !adj(P[i], P[j])) {
        mate = j;
        break;
      }
    }
    if (mate == P.size()) {
      forced.push_back(P[i]);
    } else {
      paired[i] = paired[mate] = 1;
      missingEdges.emplace_back(P[i], P[mate]);
    }
  }
}

// Stop-after-one Bron--Kerbosch for a formal branch B=(M,Q).  The initial X
// is the part that the legacy implementation omitted: common neighbors of M
// that lie outside Q.  Carrying X through the recursion makes every returned
// leaf globally maximal, not merely maximal inside M union Q.
bool ReorderSib::findOnePure(const vector<ui> &M, const vector<ui> &Q,
                             vector<ui> &found) {
  found.clear();
  if (M.empty())
    return false;

  vector<ui> common = adjList[M[0]];
  vector<ui> scratch;
  for (ui i = 1; i < (ui)M.size() && !common.empty(); i++) {
    intersectInto(scratch, common, adjList[M[i]]);
    common.swap(scratch);
  }

  vector<ui> X;
  setDiffInto(X, common, Q);
  vector<ui> R = M;
  return findOnePureRecursive(R, Q, std::move(X), found);
}

bool ReorderSib::findOnePureRecursive(vector<ui> &R, vector<ui> P, vector<ui> X,
                                      vector<ui> &found) {
  incrementSearchStateOrThrow(checksCount);
  if (R.size() + P.size() < minCliqueSize)
    return false;

  // These structural terminals are intentionally unconditional in the Pure
  // PXR lane. They do not depend on the Hybrid graph-level portfolio.
  if (P.empty()) {
    if (X.empty()) {
      found = R;
      sort(found.begin(), found.end());
      return true;
    }
    return false;
  }

  // A zero/one-candidate child can be decided without another BK level.
  if (P.size() == 1) {
    const ui extension = P.front();
    for (ui x : X)
      if (adj(x, extension))
        return false;
    found = R;
    found.push_back(extension);
    sort(found.begin(), found.end());
    return found.size() >= minCliqueSize;
  }

  ui pivot = P.front();
  ui minPScore = static_cast<ui>(P.size());
  ui universalP = numeric_limits<ui>::max();
  bool xUniversal = false;
  scanPurePXRState(P, X, pivot, minPScore, universalP, xUniversal);
  const ui pSize = static_cast<ui>(P.size());

  // An excluded vertex covering all of P makes every continuation nonmaximal.
  if (xUniversal)
    return false;

  // P is complete. With no X-universal blocker, R union P is the sole
  // maximal continuation even when X itself is nonempty.
  if (minPScore + 1 == pSize) {
    found = R;
    found.insert(found.end(), P.begin(), P.end());
    sort(found.begin(), found.end());
    return found.size() >= minCliqueSize;
  }

  // The complement of P is a matching: isolated complement vertices are
  // forced and one endpoint from every missing edge gives a maximal clique.
  if (X.empty() && minPScore + 2 >= pSize) {
    vector<ui> forced;
    vector<pair<ui, ui>> missingEdges;
    pureMatchingParts(P, forced, missingEdges);
    found = R;
    found.insert(found.end(), forced.begin(), forced.end());
    for (const auto &edge : missingEdges)
      found.push_back(edge.first);
    sort(found.begin(), found.end());
    return found.size() >= minCliqueSize;
  }

  // If every P vertex misses at most two P-neighbors, the complement consists
  // of paths and cycles. Reuse the exact 3-plex DP to obtain one witness.
  if (X.empty() && minPScore + 3 >= pSize) {
    FastPlex3Result plex = solveFastPlex3Subtree(
        adjSet, P, static_cast<ui>(R.size()), &R, minCliqueSize, nullptr);
    if (plex.handled) {
      found = std::move(plex.witness);
      sort(found.begin(), found.end());
      return plex.found;
    }
  }

  // A P-universal vertex must be present in every maximal continuation, so
  // force it into R and make a single recursive call.
  if (universalP != numeric_limits<ui>::max()) {
    vector<ui> childP;
    vector<ui> childX;
    intersectInto(childP, P, adjList[universalP]);
    intersectInto(childX, X, adjList[universalP]);
    R.push_back(universalP);
    const bool result =
        findOnePureRecursive(R, std::move(childP), std::move(childX), found);
    R.pop_back();
    return result;
  }

  vector<ui> branchRoots;
  setDiffInto(branchRoots, P, adjList[pivot]);
  for (ui v : branchRoots) {
    vector<ui> childP;
    vector<ui> childX;
    intersectInto(childP, P, adjList[v]);
    intersectInto(childX, X, adjList[v]);
    R.push_back(v);
    if (findOnePureRecursive(R, std::move(childP), std::move(childX), found)) {
      R.pop_back();
      return true;
    }
    R.pop_back();

    auto pIt = lower_bound(P.begin(), P.end(), v);
    if (pIt != P.end() && *pIt == v)
      P.erase(pIt);
    X.insert(lower_bound(X.begin(), X.end(), v), v);
  }
  return false;
}

// Exhaustive pivot Bron--Kerbosch for a formal Pure branch B=(M,Q).  Unlike
// findOnePure this consumes the whole branch in one pass.  It is used only
// when a fixed-size hitting-set mask remains over capacity after exact
// preprocessing.  The global output guard safely absorbs overlap with cliques
// already emitted by sibling-effect branches.
void ReorderSib::enumerateAllPureBranch(const vector<ui> &M,
                                        const vector<ui> &Q) {
  if (M.empty())
    return;

  vector<ui> common = adjList[M[0]];
  vector<ui> scratch;
  for (ui i = 1; i < (ui)M.size() && !common.empty(); i++) {
    intersectInto(scratch, common, adjList[M[i]]);
    common.swap(scratch);
  }

  vector<ui> X;
  setDiffInto(X, common, Q);
  vector<ui> R = M;
  enumerateAllPureBranchRecursive(R, Q, std::move(X));
}

void ReorderSib::enumerateAllPureBranchRecursive(vector<ui> &R, vector<ui> P,
                                                 vector<ui> X) {
  incrementSearchStateOrThrow(checksCount);
  if (R.size() + P.size() < minCliqueSize)
    return;

  // The same always-on terminals used by witness search also consume an
  // over-capacity Pure branch without descending through ordinary Pivot-BK.
  if (P.empty()) {
    if (X.empty())
      recordPureClique(R);
    return;
  }

  if (P.size() == 1) {
    const ui extension = P.front();
    for (ui x : X)
      if (adj(x, extension))
        return;
    vector<ui> clique = R;
    clique.push_back(extension);
    if (clique.size() >= minCliqueSize)
      recordPureClique(std::move(clique));
    return;
  }

  ui pivot = P.front();
  ui minPScore = static_cast<ui>(P.size());
  ui universalP = numeric_limits<ui>::max();
  bool xUniversal = false;
  scanPurePXRState(P, X, pivot, minPScore, universalP, xUniversal);
  const ui pSize = static_cast<ui>(P.size());

  if (xUniversal)
    return;

  if (minPScore + 1 == pSize) {
    vector<ui> clique = R;
    clique.insert(clique.end(), P.begin(), P.end());
    if (clique.size() >= minCliqueSize)
      recordPureClique(std::move(clique));
    return;
  }

  if (X.empty() && minPScore + 2 >= pSize) {
    vector<ui> forced;
    vector<pair<ui, ui>> missingEdges;
    pureMatchingParts(P, forced, missingEdges);
    vector<ui> clique = R;
    clique.insert(clique.end(), forced.begin(), forced.end());
    function<void(size_t)> materialize = [&](size_t at) {
      if (at == missingEdges.size()) {
        if (clique.size() >= minCliqueSize)
          recordPureClique(clique);
        return;
      }
      clique.push_back(missingEdges[at].first);
      materialize(at + 1);
      clique.back() = missingEdges[at].second;
      materialize(at + 1);
      clique.pop_back();
    };
    materialize(0);
    return;
  }

  if (X.empty() && minPScore + 3 >= pSize) {
    FastCliqueSink sink = [&](const vector<ui> &clique) {
      recordPureClique(clique);
    };
    FastPlex3Result plex = solveFastPlex3Subtree(
        adjSet, P, static_cast<ui>(R.size()), &R, minCliqueSize, &sink);
    if (plex.handled)
      return;
  }

  if (universalP != numeric_limits<ui>::max()) {
    vector<ui> childP;
    vector<ui> childX;
    intersectInto(childP, P, adjList[universalP]);
    intersectInto(childX, X, adjList[universalP]);
    R.push_back(universalP);
    enumerateAllPureBranchRecursive(R, std::move(childP), std::move(childX));
    R.pop_back();
    return;
  }

  vector<ui> branchRoots;
  setDiffInto(branchRoots, P, adjList[pivot]);
  for (ui v : branchRoots) {
    vector<ui> childP;
    vector<ui> childX;
    intersectInto(childP, P, adjList[v]);
    intersectInto(childX, X, adjList[v]);
    R.push_back(v);
    enumerateAllPureBranchRecursive(R, std::move(childP), std::move(childX));
    R.pop_back();

    auto pIt = lower_bound(P.begin(), P.end(), v);
    if (pIt != P.end() && *pIt == v)
      P.erase(pIt);
    X.insert(lower_bound(X.begin(), X.end(), v), v);
  }
}

static string encodeClique(const vector<ui> &C) {
  ScopedTimer _t(rsp.encode_ms, rsp.encode_n);
  return string(reinterpret_cast<const char *>(C.data()),
                C.size() * sizeof(ui));
}

bool ReorderSib::recordPureClique(vector<ui> C) {
  sort(C.begin(), C.end());
  const string key = encodeClique(C);
  if (emittedCliqueKeys.find(key) != emittedCliqueKeys.end()) {
    addCliqueCountOrThrow(dupBlocked, 1);
    return false;
  }

  if (allCliques.size() > numeric_limits<ui>::max())
    throw overflow_error("materialized clique index exceeds uint32_t");
  for (ui v : C) {
    if (cliqueCountByVertex[v] == numeric_limits<ull>::max())
      throw overflow_error("per-vertex clique count exceeds uint64_t");
  }

  const ui cliqueId = static_cast<ui>(allCliques.size());
  addCliqueCountOrThrow(cliqueCount, 1);
  emittedCliqueKeys.insert(key);
  allCliques.push_back(std::move(C));
  maxCliqueSize = max(maxCliqueSize, allCliques.back().size());
  for (ui v : allCliques.back()) {
    if (cliquesByVertexByLevel[v].empty())
      cliquesByVertexByLevel[v].resize(1);
    cliquesByVertexByLevel[v][0].push_back(cliqueId);
    addCliqueCountOrThrow(cliqueCountByVertex[v], 1);
  }
  return true;
}

vector<vector<ui>> ReorderSib::getCliques() const {
  vector<vector<ui>> restored = allCliques;
  for (vector<ui> &clique : restored) {
    for (ui &vertex : clique)
      vertex = internalToOriginal[vertex];
    sort(clique.begin(), clique.end());
  }
  return restored;
}

void ReorderSib::findAllMaximalCliquesPure() {
  rsp.reset();
  cliqueCount = 0;
  dupBlocked = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  solverBudgetFallbacks = 0;
  allCliques.clear();
  emittedCliqueKeys.clear();
  cliquesByVertexByLevel.assign(n, {});
  fill(cliqueCountByVertex.begin(), cliqueCountByVertex.end(), 0);

  vector<PureBranch> worklist;
  worklist.reserve(n);
  // Reverse insertion makes the stack visit canonical roots from low to high.
  for (ui next = n; next > 0; next--) {
    const ui v = next - 1;
    worklist.push_back({{v}, adjList2[v]});
  }

  auto t0 = chrono::high_resolution_clock::now();
  while (!worklist.empty()) {
    PureBranch branch = std::move(worklist.back());
    worklist.pop_back();

    if (branch.mustin.size() + branch.expandTo.size() < minCliqueSize)
      continue;

    vector<ui> covers = collectAllCoveringCliques(branch.mustin);
    if (covers.empty()) {
      vector<ui> found;
      if (findOnePure(branch.mustin, branch.expandTo, found)) {
        recordPureClique(std::move(found));
        // The unchanged branch retains every unseen target.  On its next pop,
        // the clique just recorded necessarily covers its must-in set.
        worklist.push_back(std::move(branch));
      }
      continue;
    }

    bool usePivotFallback = false;
    vector<vector<ui>> seeds =
        generateExactSiblingSets(branch.expandTo, covers, &usePivotFallback);
    if (usePivotFallback) {
      enumerateAllPureBranch(branch.mustin, branch.expandTo);
      continue;
    }
    // Reverse insertion preserves the solver's deterministic seed order under
    // the LIFO worklist; correctness does not depend on this order.
    for (auto it = seeds.rbegin(); it != seeds.rend(); ++it) {
      vector<ui> nextM = unionSet(branch.mustin, *it);
      vector<ui> nextQ = commonExpand(branch.expandTo, *it);
      worklist.push_back({std::move(nextM), std::move(nextQ)});
    }
  }
  auto t1 = chrono::high_resolution_clock::now();
  const double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << fixed << setprecision(3) << "PureReorderSib: cliques=" << cliqueCount
       << "  dups=" << dupBlocked << "  maxSize=" << maxCliqueSize
       << "  minSize=" << minCliqueSize << "  checks=" << checksCount
       << "  budgetFallbacks=" << solverBudgetFallbacks << "  time=" << ms
       << " ms" << endl;
#if PROFILING
  rsp.print(ms);
#endif
}
