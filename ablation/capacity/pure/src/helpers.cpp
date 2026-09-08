#include "../inc/helpers.h"
#include "../inc/config.h"
#include "../inc/fast_plex3.h"
#if !defined(PURE_LEAN_BENCHMARK)
#include <chrono>
#include <iomanip>
#endif
#include <functional>
#include <numeric>
#include <type_traits>

namespace {
using pure_config::kEt1Enabled;
using pure_config::kEt2Enabled;
using pure_config::kEt3Enabled;
using pure_config::kAdjHashThreshold;
using pure_config::kDiagnosticsEnabled;
using pure_config::kHitsetCapacity;
using pure_config::kHitsetDynamic;
using pure_config::kPruneAntichain;
using pure_config::kPruneFailFirst;
using pure_config::kPruneNormalization;
using pure_config::kPruneSubsumption;
using pure_config::kPruneUnit;
using pure_config::kPruneUsefulness;
using pure_config::kPruneZeroCoverage;

constexpr size_t kBinaryIntersectionRatio = 16;

vector<vector<ui>>
minimalByInclusion(vector<vector<ui>> solutions) {
  for (vector<ui> &solution : solutions)
    sort(solution.begin(), solution.end());
  sort(solutions.begin(), solutions.end());
  solutions.erase(unique(solutions.begin(), solutions.end()),
                  solutions.end());
  vector<vector<ui>> minimal;
  for (size_t i = 0; i < solutions.size(); ++i) {
    bool hasSmallerCover = false;
    for (size_t j = 0; j < solutions.size(); ++j) {
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
} // namespace

// Returns peelSeq index of verticies by core value
// peelSeq[0] = highest-core vertex, peelSeq[n-1] = lowest.
static vector<ui> computePeelSeq(const Graph &g, ui *degeneracy = nullptr) {
  ui n = g.n;
  if (degeneracy != nullptr)
    *degeneracy = 0;
  if (n == 0)
    return {};

  vector<ui> deg(g.degree.begin(), g.degree.end());
  ui maxDeg = *max_element(deg.begin(), deg.end());
  ui maxCore = 0;

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
    maxCore = max(maxCore, deg[v]);
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
  if (degeneracy != nullptr)
    *degeneracy = maxCore;
  return peelSeq;
}

// Build the reordered graph directly as sorted CSR and remember the local
// offset at which each row's higher-index neighbors begin.
static void buildAdjacencyCsr(const Graph &g, const vector<ui> &perm,
                              vector<ui> &adjVertices,
                              vector<size_t> &adjOffsets,
                              vector<ui> &firstForwardNeighbor) {
  const ui n = g.n;
  adjOffsets.assign(static_cast<size_t>(n) + 1, 0);
  firstForwardNeighbor.assign(n, 0);
  for (ui original = 0; original < n; ++original) {
    const ui reordered = perm[original];
    adjOffsets[static_cast<size_t>(reordered) + 1] =
        static_cast<size_t>(g.offset[original + 1] - g.offset[original]);
  }
  for (ui reordered = 0; reordered < n; ++reordered)
    adjOffsets[static_cast<size_t>(reordered) + 1] += adjOffsets[reordered];

  adjVertices.resize(adjOffsets[n]);
  for (ui original = 0; original < n; ++original) {
    const ui reordered = perm[original];
    auto rowBegin = adjVertices.begin() + adjOffsets[reordered];
    size_t at = 0;
    for (ui j = g.offset[original]; j < g.offset[original + 1]; ++j)
      rowBegin[at++] = perm[g.neighbors[j]];
    auto rowEnd = adjVertices.begin() + adjOffsets[reordered + 1];
    sort(rowBegin, rowEnd);
    firstForwardNeighbor[reordered] = static_cast<ui>(
        upper_bound(rowBegin, rowEnd, reordered) - rowBegin);
  }
}

// ReorderSib Implementation
ReorderSib::ReorderSib(Graph &g, SibMethod method, ui minCliqueSize)
    : method(method), minCliqueSize(max<ui>(1, minCliqueSize)) {
  n = g.n;
  cliqueCount = 0;
  dupBlocked = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  findOnePxrStates = 0;
  fullPxrStates = 0;
  et1EnumeratedStates = 0;
  et2EnumeratedStates = 0;
  et3EnumeratedStates = 0;
  solverWorkBudget = 0;
  solverWorkBudgetEnabled = false;
  solverBudgetFallbacks = 0;
  solverCertifiedBudgetFallbacks = 0;
  solverCapacityFallbacks = 0;
  nontrivialSeedSolverCalls = 0;
  maximumSeedConstraints = 0;
  findOne2PlexTerminals = 0;
  findOne3PlexTerminals = 0;
  fullPxr2PlexTerminals = 0;
  fullPxr3PlexTerminals = 0;
  worklistPushes = 0;
  worklistPops = 0;
  maximumWorklistSize = 0;
  minSizePrunedBranches = 0;
  coverLookupCalls = 0;
  findOneCalls = 0;
  findOneSuccesses = 0;
  findOneCliqueSizeTotal = 0;
  maximumFindOneCliqueSize = 0;
  seedSolverCalls = 0;
  generatedBranches = 0;
  fullPxrFallbackBranches = 0;
  smallQFullPxrBranches = 0;
  poppedQSizeBuckets.fill(0);
  smallQFullPxrThreshold = 4;

  cliqueIdsByVertex.resize(n);

  vector<ui> perm(n);

  vector<ui> peelSeq = computePeelSeq(g);
  for (ui i = 0; i < n; i++)
    perm[peelSeq[n - 1 - i]] = i; // low-core → low index

  internalToOriginal.resize(n);
  for (ui original = 0; original < n; original++)
    internalToOriginal[perm[original]] = original;

  buildAdjacencyCsr(g, perm, adjVertices, adjOffsets,
                    firstForwardNeighbor);

  adjHash.resize(n);
  for (ui u = 0; u < n; u++) {
    const AdjacencyRow row = adjacentVertices(u);
    if (row.size() < kAdjHashThreshold)
      continue;
    auto hash = make_unique<unordered_set<ui>>();
    hash->reserve(row.size());
    hash->insert(row.begin(), row.end());
    adjHash[u] = std::move(hash);
  }

  eIndex.assign(n, 0);
  eIndexStamp.assign(n, 0);
  eIndexToken = 0;
}

void ReorderSib::intersectInto(vector<ui> &out, const vector<ui> &A,
                               const vector<ui> &B) {
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    auto search = B.begin();
    for (ui v : A) {
      search = lower_bound(search, B.end(), v);
      if (search == B.end())
        break;
      if (*search == v) {
        out.push_back(v);
        ++search;
      }
    }
    return;
  }
  if (B.size() * kBinaryIntersectionRatio < A.size()) {
    auto search = A.begin();
    for (ui v : B) {
      search = lower_bound(search, A.end(), v);
      if (search == A.end())
        break;
      if (*search == v) {
        out.push_back(v);
        ++search;
      }
    }
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

void ReorderSib::intersectInto(vector<ui> &out, const vector<ui> &A,
                               AdjacencyRow B) {
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    auto search = B.begin();
    for (ui v : A) {
      search = lower_bound(search, B.end(), v);
      if (search == B.end())
        break;
      if (*search == v) {
        out.push_back(v);
        ++search;
      }
    }
    return;
  }
  if (B.size() * kBinaryIntersectionRatio < A.size()) {
    auto search = A.begin();
    for (ui v : B) {
      search = lower_bound(search, A.end(), v);
      if (search == A.end())
        break;
      if (*search == v) {
        out.push_back(v);
        ++search;
      }
    }
    return;
  }
  size_t i = 0, j = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] == B[j]) {
      out.push_back(A[i]);
      ++i;
      ++j;
    } else if (A[i] < B[j]) {
      ++i;
    } else {
      ++j;
    }
  }
}

void ReorderSib::intersectExcludingInto(vector<ui> &out, const vector<ui> &A,
                                        const vector<ui> &B,
                                        const vector<ui> &exclude) {
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    for (ui v : A) {
      if (binary_search(B.begin(), B.end(), v) &&
          !binary_search(exclude.begin(), exclude.end(), v))
        out.push_back(v);
    }
    return;
  }
  if (B.size() * kBinaryIntersectionRatio < A.size()) {
    for (ui v : B) {
      if (binary_search(A.begin(), A.end(), v) &&
          !binary_search(exclude.begin(), exclude.end(), v))
        out.push_back(v);
    }
    return;
  }
  ui i = 0, j = 0, k = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] == B[j]) {
      while (k < exclude.size() && exclude[k] < A[i])
        k++;
      if (k == exclude.size() || exclude[k] != A[i])
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

void ReorderSib::intersectExcludingInto(vector<ui> &out,
                                        const vector<ui> &A, AdjacencyRow B,
                                        const vector<ui> &exclude) {
  out.clear();
  const size_t need = min(A.size(), B.size());
  if (out.capacity() < need)
    out.reserve(need);
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    for (ui v : A) {
      if (binary_search(B.begin(), B.end(), v) &&
          !binary_search(exclude.begin(), exclude.end(), v))
        out.push_back(v);
    }
    return;
  }
  if (B.size() * kBinaryIntersectionRatio < A.size()) {
    for (ui v : B) {
      if (binary_search(A.begin(), A.end(), v) &&
          !binary_search(exclude.begin(), exclude.end(), v))
        out.push_back(v);
    }
    return;
  }
  size_t i = 0, j = 0, k = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] == B[j]) {
      while (k < exclude.size() && exclude[k] < A[i])
        ++k;
      if (k == exclude.size() || exclude[k] != A[i])
        out.push_back(A[i]);
      ++i;
      ++j;
    } else if (A[i] < B[j]) {
      ++i;
    } else {
      ++j;
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

vector<ui> ReorderSib::setDiffStoredClique(const vector<ui> &A,
                                           ui cliqueId) const {
  vector<ui> result;
  result.reserve(A.size());
  size_t i = 0;
  size_t j = cliqueOffsets[cliqueId];
  const size_t end = cliqueOffsets[cliqueId + 1];
  while (i < A.size()) {
    if (j == end || A[i] < cliqueVertices[j]) {
      result.push_back(A[i++]);
    } else if (A[i] == cliqueVertices[j]) {
      ++i;
      ++j;
    } else {
      ++j;
    }
  }
  return result;
}

bool ReorderSib::storedCliqueContains(
    ui cliqueId, const vector<ui> &subset) const {
  const auto begin = cliqueVertices.begin() + cliqueOffsets[cliqueId];
  const auto end = cliqueVertices.begin() + cliqueOffsets[cliqueId + 1];
  return includes(begin, end, subset.begin(), subset.end());
}

bool ReorderSib::storedCliqueEquals(ui cliqueId,
                                    const vector<ui> &clique) const {
  if (storedCliqueSize(cliqueId) != clique.size())
    return false;
  const auto begin = cliqueVertices.begin() + cliqueOffsets[cliqueId];
  const auto end = cliqueVertices.begin() + cliqueOffsets[cliqueId + 1];
  return equal(begin, end, clique.begin());
}

void ReorderSib::setDiffInto(vector<ui> &out, const vector<ui> &A,
                             const vector<ui> &B) {
  out.clear();
  if (out.capacity() < A.size())
    out.reserve(A.size());
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    auto search = B.begin();
    const auto end = B.end();
    for (ui v : A) {
      search = lower_bound(search, end, v);
      if (search == end || *search != v)
        out.push_back(v);
      else
        ++search;
    }
    return;
  }
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

void ReorderSib::setDiffInto(vector<ui> &out, const vector<ui> &A,
                             AdjacencyRow B) {
  out.clear();
  if (out.capacity() < A.size())
    out.reserve(A.size());
  if (A.size() * kBinaryIntersectionRatio < B.size()) {
    auto search = B.begin();
    const auto end = B.end();
    for (ui v : A) {
      search = lower_bound(search, end, v);
      if (search == end || *search != v)
        out.push_back(v);
      else
        ++search;
    }
    return;
  }
  size_t i = 0, j = 0;
  while (i < A.size()) {
    if (j == B.size() || A[i] < B[j]) {
      out.push_back(A[i]);
      ++i;
    } else if (A[i] == B[j]) {
      ++i;
      ++j;
    } else {
      ++j;
    }
  }
}

void ReorderSib::unionSetInto(vector<ui> &out, const vector<ui> &A,
                              const vector<ui> &B) {
  out.clear();
  if (out.capacity() < A.size() + B.size())
    out.reserve(A.size() + B.size());
  ui i = 0, j = 0;
  while (i < A.size() && j < B.size()) {
    if (A[i] < B[j]) {
      out.push_back(A[i]);
      i++;
    } else if (A[i] > B[j]) {
      out.push_back(B[j]);
      j++;
    } else {
      out.push_back(A[i]);
      i++;
      j++;
    }
  }
  while (i < A.size())
    out.push_back(A[i++]);
  while (j < B.size())
    out.push_back(B[j++]);
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
void ReorderSib::commonExpandInto(vector<ui> &out, const vector<ui> &E,
                                  const vector<ui> &S) {
  if (S.empty()) {
    out.assign(E.begin(), E.end());
    return;
  }
  if (S.size() == 1) {
    intersectInto(out, E, adjacentVertices(S[0]));
    return;
  }

  commonExpandOrder.assign(S.begin(), S.end());
  sort(commonExpandOrder.begin(), commonExpandOrder.end(), [&](ui a, ui b) {
    return adjacentVertices(a).size() < adjacentVertices(b).size();
  });

  intersectExcludingInto(out, E, adjacentVertices(commonExpandOrder[0]), S);

  for (ui idx = 1;
       idx < (ui)commonExpandOrder.size() && !out.empty(); idx++) {
    intersectInto(commonExpandScratch, out,
                  adjacentVertices(commonExpandOrder[idx]));
    out.swap(commonExpandScratch);
  }
}

// Proof-faithful cover lookup for Pure-ReorderSib. Every emitted clique
// containing M must be visible before FindOne may run, so no recent-level or
// weak-cover filters are applied here.
vector<ui> ReorderSib::collectAllCoveringCliques(const vector<ui> &M) {
  vector<ui> result;
  if (M.empty())
    return result;

  ui seed = M[0];
  for (ui v : M)
    if (cliqueIdsByVertex[v].size() < cliqueIdsByVertex[seed].size())
      seed = v;

  const auto &posting = cliqueIdsByVertex[seed];
  result.reserve(posting.size());
  for (ui cliqueId : posting) {
    if (storedCliqueContains(cliqueId, M))
      result.push_back(cliqueId);
  }
  return result;
}

// Each covering clique C contributes the constraint "pick something from E
// that is outside C", which is exactly E \ C.
vector<vector<ui>> ReorderSib::buildHitSets(const vector<ui> &E,
                                            const vector<ui> &cliqueIds,
                                            ui maxHitSets) {
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
      return storedCliqueSize(a) > storedCliqueSize(b);
    });
    sorted.resize(maxHitSets);
    ids = &sorted;
  }

  vector<vector<ui>> hitSets;
  hitSets.reserve(ids->size());
  for (ui cId : *ids)
    hitSets.push_back(setDiffStoredClique(E, cId));
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

// Exact sibling-seed generation used by the theorem-aligned Pure lane. It
// applies every covering-clique constraint.
vector<vector<ui>> ReorderSib::generateExactSiblingSets(
    const vector<ui> &E, const vector<ui> &coveringCliqueIds,
    bool *usePivotFallback) {
  if (usePivotFallback != nullptr)
    *usePivotFallback = false;
  if (coveringCliqueIds.empty())
    return {};

  // Keep the direct compact-arena fast path identical across every capacity
  // variant. Larger inputs go through the shared representation-independent
  // preprocessing before the configured fixed/dynamic mask is selected.
  if (method == SibMethod::OPTIMIZED &&
      coveringCliqueIds.size() <= 64)
    return efficientHittingSetDirect(E, coveringCliqueIds,
                                     usePivotFallback);

  vector<vector<ui>> hitSets = buildHitSets(E, coveringCliqueIds);
  for (const vector<ui> &hitSet : hitSets)
    if (hitSet.empty())
      return {};

  if (coveringCliqueIds.size() == 1)
    return singletonBranches(hitSets[0]);
  if (method == SibMethod::BACKTRACKING)
    return backtrackingBranchBound(E, hitSets);
  return efficientHittingSet(E, std::move(hitSets), usePivotFallback);
}

vector<vector<ui>>
ReorderSib::backtrackingBranchBound(const vector<ui> &E,
                                    const vector<vector<ui>> &hitSets) {
  vector<vector<ui>> solutions;
  vector<ui> current;

  // DFS over clique-compatible subsets of E. Once the current set already hits
  // every constraint, maintain the minimal solution family online and stop
  // descending that branch.
  function<void(ui)> dfs = [&](ui start) {
    if (hitsAll(current, hitSets)) {
      // current is already sorted because DFS only appends E[i] in increasing
      // index order, and E itself is sorted.
      for (const auto &s : solutions)
        if (includes(current.begin(), current.end(), s.begin(), s.end()))
          return;

      solutions.erase(remove_if(solutions.begin(), solutions.end(),
                                [&](const vector<ui> &s) {
                                  return includes(s.begin(), s.end(),
                                                  current.begin(),
                                                  current.end());
                                }),
                      solutions.end());
      solutions.push_back(current);
      return;
    }

    // If current already contains a known minimal solution, any extension is a
    // non-minimal superset and can be pruned immediately.
    for (const auto &s : solutions)
      if (includes(current.begin(), current.end(), s.begin(), s.end()))
        return;

    for (ui i = start; i < E.size(); i++) {
      bool connected = true;
      for (ui v : current) {
        const AdjacencyRow row = adjacentVertices(v);
        if (!binary_search(row.begin(), row.end(), E[i])) {
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
  return solutions;
}

// Construct the fixed-width seed-solver coverage masks straight from the
// compact clique arena.  This avoids allocating one vector for every E \\ C
// constraint on the overwhelmingly common <=128-constraint path.
vector<vector<ui>> ReorderSib::efficientHittingSetDirect(
    const vector<ui> &inputE, const vector<ui> &coveringCliqueIds,
    bool *usePivotFallback) {
  if (usePivotFallback != nullptr)
    *usePivotFallback = false;

  static constexpr ui fixedMaskWords = 1;
  using RawMask = array<ull, fixedMaskWords>;
  auto zeroRawMask = []() -> RawMask { return RawMask{}; };
  auto setRawBit = [](RawMask &m, ui bit) {
    m[bit >> 6] |= (1ULL << (bit & 63));
  };
  auto hasRawBit = [](const RawMask &m, ui bit) -> bool {
    return (m[bit >> 6] & (1ULL << (bit & 63))) != 0;
  };

  const ui rawHSize = static_cast<ui>(coveringCliqueIds.size());
  vector<RawMask> rawCoverage(inputE.size(), zeroRawMask());
  vector<ui> constraintSizes(rawHSize, 0);
  vector<ull> constraintHashes(rawHSize, 1469598103934665603ULL);

  for (ui bit = 0; bit < rawHSize; ++bit) {
    const ui cliqueId = coveringCliqueIds[bit];
    size_t cliqueAt = cliqueOffsets[cliqueId];
    const size_t cliqueEnd = cliqueOffsets[cliqueId + 1];
    for (size_t i = 0; i < inputE.size(); ++i) {
      const ui vertex = inputE[i];
      while (cliqueAt < cliqueEnd && cliqueVertices[cliqueAt] < vertex)
        ++cliqueAt;
      if (cliqueAt < cliqueEnd && cliqueVertices[cliqueAt] == vertex)
        continue;
      setRawBit(rawCoverage[i], bit);
      ++constraintSizes[bit];
      constraintHashes[bit] ^= static_cast<ull>(vertex);
      constraintHashes[bit] *= 1099511628211ULL;
    }
    if (constraintSizes[bit] == 0)
      return {};
    constraintHashes[bit] ^= constraintSizes[bit];
    constraintHashes[bit] *= 1099511628211ULL;
  }

  auto sameConstraint = [&](ui a, ui b) {
    for (const RawMask &coverage : rawCoverage)
      if (hasRawBit(coverage, a) != hasRawBit(coverage, b))
        return false;
    return true;
  };
  vector<ui> constraintOrder;
  constraintOrder.reserve(rawHSize);
  for (ui bit = 0; bit < rawHSize; ++bit) {
    bool duplicate = false;
    if (kPruneNormalization) {
      for (ui kept : constraintOrder) {
        if (constraintSizes[kept] == constraintSizes[bit] &&
            constraintHashes[kept] == constraintHashes[bit] &&
            sameConstraint(kept, bit)) {
          duplicate = true;
          break;
        }
      }
    }
    if (!duplicate)
      constraintOrder.push_back(bit);
  }
  sort(constraintOrder.begin(), constraintOrder.end(), [&](ui a, ui b) {
    if (constraintSizes[a] != constraintSizes[b])
      return constraintSizes[a] < constraintSizes[b];
    for (const RawMask &coverage : rawCoverage) {
      const bool inA = hasRawBit(coverage, a);
      const bool inB = hasRawBit(coverage, b);
      if (inA != inB)
        return inA;
    }
    return a < b;
  });

  const ui hSize = static_cast<ui>(constraintOrder.size());
  addCliqueCountOrThrow(nontrivialSeedSolverCalls, 1);
  maximumSeedConstraints = max(maximumSeedConstraints, hSize);
  auto solve = [&](auto wordCountTag) -> vector<vector<ui>> {
  static constexpr ui maskWords = decltype(wordCountTag)::value;
  using Mask = array<ull, maskWords>;
  auto zeroMask = []() -> Mask { return Mask{}; };
  auto setBit = [](Mask &m, ui bit) {
    m[bit >> 6] |= (1ULL << (bit & 63));
  };
  auto clearBit = [](Mask &m, ui bit) {
    m[bit >> 6] &= ~(1ULL << (bit & 63));
  };
  auto hasBit = [](const Mask &m, ui bit) -> bool {
    return (m[bit >> 6] & (1ULL << (bit & 63))) != 0;
  };
  auto anyMask = [](const Mask &m) -> bool {
    for (ui word = 0; word < maskWords; ++word)
      if (m[word] != 0)
        return true;
    return false;
  };
  auto orEq = [](Mask &dst, const Mask &src) {
    for (ui word = 0; word < maskWords; ++word)
      dst[word] |= src[word];
  };
  auto andNot = [&](const Mask &a, const Mask &b) -> Mask {
    Mask result = zeroMask();
    for (ui word = 0; word < maskWords; ++word)
      result[word] = a[word] & ~b[word];
    return result;
  };
  auto intersects = [](const Mask &a, const Mask &b) -> bool {
    for (ui word = 0; word < maskWords; ++word)
      if ((a[word] & b[word]) != 0)
        return true;
    return false;
  };
  auto coversAll = [](const Mask &covered, const Mask &need) -> bool {
    for (ui word = 0; word < maskWords; ++word)
      if ((covered[word] & need[word]) != need[word])
        return false;
    return true;
  };
  auto popcountMask = [](const Mask &m) -> int {
    int count = 0;
    for (ui word = 0; word < maskWords; ++word)
      count += __builtin_popcountll(m[word]);
    return count;
  };
  auto popNextBit = [](Mask &m) -> int {
    for (ui word = 0; word < maskWords; ++word) {
      if (m[word] == 0)
        continue;
      const int bit = __builtin_ctzll(m[word]);
      m[word] &= m[word] - 1;
      return static_cast<int>(word * 64 + bit);
    }
    return -1;
  };

  Mask fullMask = zeroMask();
  for (ui bit = 0; bit < hSize; ++bit)
    setBit(fullMask, bit);

  vector<ui> E;
  vector<Mask> cov;
  E.reserve(inputE.size());
  cov.reserve(inputE.size());
  for (size_t i = 0; i < inputE.size(); ++i) {
    Mask dense = zeroMask();
    for (ui bit = 0; bit < hSize; ++bit)
      if (hasRawBit(rawCoverage[i], constraintOrder[bit]))
        setBit(dense, bit);
    if (!kPruneUsefulness || anyMask(dense)) {
      E.push_back(inputE[i]);
      cov.push_back(dense);
    }
  }
  if (E.empty())
    return {};

  const ui eSize = static_cast<ui>(E.size());

  vector<ui> initCands(eSize);
  iota(initCands.begin(), initCands.end(), 0);
  sort(initCands.begin(), initCands.end(), [&](ui a, ui b) {
    const int aCoverage = popcountMask(cov[a]);
    const int bCoverage = popcountMask(cov[b]);
    if (aCoverage != bCoverage)
      return aCoverage > bCoverage;
    return E[a] < E[b];
  });

  vector<ui> forcedIdxs;
  Mask forcedCov = zeroMask();
  if (kPruneUnit) {
    vector<ui> activeCands = initCands;
    bool changed = true;
    bool conflict = false;
    while (changed && !conflict) {
      changed = false;
      Mask remaining = andNot(fullMask, forcedCov);
      if (!anyMask(remaining))
        break;
      while (!conflict) {
        const int h = popNextBit(remaining);
        if (h < 0)
          break;
        ui sole = eSize;
        int count = 0;
        for (ui candidate : activeCands) {
          if (!hasBit(cov[candidate], static_cast<ui>(h)))
            continue;
          sole = candidate;
          if (++count > 1)
            break;
        }
        if (count == 0) {
          conflict = true;
          break;
        }
        if (count != 1)
          continue;
        for (ui forced : forcedIdxs) {
          if (!adj(E[forced], E[sole])) {
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
        for (ui candidate : activeCands)
          if (candidate != sole && adj(E[sole], E[candidate]))
            next.push_back(candidate);
        activeCands = std::move(next);
        changed = true;
        break;
      }
    }
    if (conflict)
      return {};
    initCands = std::move(activeCands);
  }

  if (kPruneSubsumption) {
    for (ui h = 0; h < hSize; ++h) {
      if (!hasBit(fullMask, h) || hasBit(forcedCov, h))
        continue;
      for (ui g = 0; g < hSize; ++g) {
        if (g == h || !hasBit(fullMask, g) || hasBit(forcedCov, g))
          continue;
        bool gSubsumesH = true;
        for (ui candidate : initCands) {
          if (hasBit(cov[candidate], g) &&
              !hasBit(cov[candidate], h)) {
            gSubsumesH = false;
            break;
          }
        }
        if (gSubsumesH) {
          clearBit(fullMask, h);
          break;
        }
      }
    }
  }

  if (coversAll(forcedCov, fullMask)) {
    if (forcedIdxs.empty())
      return {};
    vector<ui> solution;
    solution.reserve(forcedIdxs.size());
    for (ui idx : forcedIdxs)
      solution.push_back(E[idx]);
    sort(solution.begin(), solution.end());
    return {std::move(solution)};
  }

  // Before paying O(|E|^2) to materialize compatibility, prove whether the
  // root alone must exceed the configured DFS pair-check budget. Every root
  // candidate that adds coverage scans every numerically later such candidate,
  // so C(k, 2) is a strict lower bound on unavoidable root work. Verify root
  // feasibility first because fail-first can otherwise reject without work.
  if (solverWorkBudgetEnabled && usePivotFallback != nullptr) {
    const ull candidateCount = static_cast<ull>(initCands.size());
    if (candidateCount > 1 &&
        candidateCount * (candidateCount - 1) / 2 > solverWorkBudget) {
      const Mask uncovered = andNot(fullMask, forcedCov);
      Mask available = forcedCov;
      ull rootCandidates = 0;
      for (ui candidate : initCands) {
        orEq(available, cov[candidate]);
        if (!kPruneZeroCoverage || intersects(cov[candidate], uncovered))
          ++rootCandidates;
      }
      if (!kPruneFailFirst || coversAll(available, fullMask)) {
        const ull unavoidableRootWork =
            rootCandidates * (rootCandidates - 1) / 2;
        if (unavoidableRootWork > solverWorkBudget) {
          *usePivotFallback = true;
          if constexpr (kDiagnosticsEnabled) {
            addCliqueCountOrThrow(solverBudgetFallbacks, 1);
            addCliqueCountOrThrow(solverCertifiedBudgetFallbacks, 1);
          }
          return {};
        }
      } else {
        return {};
      }
    }
  }

  vector<char> compat(static_cast<size_t>(eSize) * eSize, 0);
  for (ui i = 0; i < eSize; ++i) {
    const AdjacencyRow row = adjacentVertices(E[i]);
    for (ui j = i + 1; j < eSize; ++j)
      if (binary_search(row.begin(), row.end(), E[j]))
        compat[static_cast<size_t>(i) * eSize + j] =
            compat[static_cast<size_t>(j) * eSize + i] = 1;
  }

  vector<vector<ui>> solutions;
  vector<ui> current;
  vector<vector<ui>> candidateScratch(eSize + 1);
  ull solverWork = 0;
  bool budgetExceeded = false;
  function<void(const vector<ui> &, Mask, ui)> dfs =
      [&](const vector<ui> &candidates, Mask covered, ui depth) {
        if (budgetExceeded)
          return;
        if (coversAll(covered, fullMask)) {
          if (kPruneAntichain) {
            for (const vector<ui> &solution : solutions)
              if (includes(current.begin(), current.end(), solution.begin(),
                           solution.end()))
                return;
            solutions.erase(
                remove_if(solutions.begin(), solutions.end(),
                          [&](const vector<ui> &solution) {
                            return includes(solution.begin(), solution.end(),
                                            current.begin(), current.end());
                          }),
                solutions.end());
          }
          solutions.push_back(current);
          return;
        }

        const Mask uncovered = andNot(fullMask, covered);
        if (kPruneAntichain)
          for (const vector<ui> &solution : solutions)
            if (includes(current.begin(), current.end(), solution.begin(),
                         solution.end()))
              return;

        if (kPruneFailFirst) {
          Mask remaining = uncovered;
          while (true) {
            const int h = popNextBit(remaining);
            if (h < 0)
              break;
            bool found = false;
            for (ui candidate : candidates)
              if (hasBit(cov[candidate], static_cast<ui>(h))) {
                found = true;
                break;
              }
            if (!found)
              return;
          }
        }

        for (ui candidate : candidates) {
          if (kPruneZeroCoverage &&
              !intersects(cov[candidate], uncovered))
            continue;
          vector<ui> &next = candidateScratch[depth];
          next.clear();
          if (next.capacity() < candidates.size())
            next.reserve(candidates.size());
          for (ui later : candidates) {
            if (later <= candidate)
              continue;
            if (solverWorkBudgetEnabled && usePivotFallback != nullptr &&
                solverWork >= solverWorkBudget) {
              budgetExceeded = true;
              break;
            }
            ++solverWork;
            if (compat[static_cast<size_t>(candidate) * eSize + later])
              next.push_back(later);
          }
          if (budgetExceeded)
            return;
          current.push_back(candidate);
          Mask nextCovered = covered;
          orEq(nextCovered, cov[candidate]);
          dfs(next, nextCovered, depth + 1);
          current.pop_back();
          if (budgetExceeded)
            return;
        }
      };

  dfs(initCands, forcedCov, 1);
  if (budgetExceeded) {
    *usePivotFallback = true;
    if constexpr (kDiagnosticsEnabled)
      addCliqueCountOrThrow(solverBudgetFallbacks, 1);
    return {};
  }
  if (!kPruneAntichain)
    solutions = minimalByInclusion(std::move(solutions));

  vector<vector<ui>> result;
  result.reserve(solutions.size());
  for (const vector<ui> &solution : solutions) {
    vector<ui> vertices;
    vertices.reserve(solution.size() + forcedIdxs.size());
    for (ui idx : solution)
      vertices.push_back(E[idx]);
    for (ui idx : forcedIdxs)
      vertices.push_back(E[idx]);
    sort(vertices.begin(), vertices.end());
    result.push_back(std::move(vertices));
  }
  return result;
  };

  if (hSize <= 64)
    return solve(integral_constant<ui, 1>{});
  return solve(integral_constant<ui, fixedMaskWords>{});
}

// Optimized exact solver for all minimal clique-constrained hitting sets.
//
// Improvements over backtrackingBranchBound:
//   1. Bitmask coverage  — done-check and update are O(1) bitwise ops.
//   2. Incremental candidate list — compat-filtered frontier passed down,
//      no per-step binary_search into reordered adjacency.
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
                                vector<vector<ui>> hitSets,
                                bool *usePivotFallback) {
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

  vector<ui> preForced;
  bool infeasible = false;
  vector<size_t> constraintOffsets(inputE.size() + 2, 0);
  vector<size_t> constraintNext;
  vector<vector<ui>> constraintSortScratch;

  auto normalizeAndSubsume = [&]() {
    // buildHitSets creates sorted unique subsequences of inputE, and every
    // filtering round below preserves that invariant.
    fill(constraintOffsets.begin(), constraintOffsets.end(), 0);
    for (const vector<ui> &H : hitSets)
      ++constraintOffsets[H.size() + 1];
    partial_sum(constraintOffsets.begin(), constraintOffsets.end(),
                constraintOffsets.begin());
    constraintNext.assign(constraintOffsets.begin(),
                          constraintOffsets.end() - 1);
    constraintSortScratch.clear();
    constraintSortScratch.resize(hitSets.size());
    for (vector<ui> &H : hitSets)
      constraintSortScratch[constraintNext[H.size()]++] = std::move(H);
    for (size_t size = 0; size + 1 < constraintOffsets.size(); ++size) {
      const size_t begin = constraintOffsets[size];
      const size_t end = constraintOffsets[size + 1];
      if (end - begin > 1)
        sort(constraintSortScratch.begin() + begin,
             constraintSortScratch.begin() + end);
    }
    hitSets.swap(constraintSortScratch);
    if (kPruneNormalization)
      hitSets.erase(unique(hitSets.begin(), hitSets.end()), hitSets.end());
    if (!kPruneSubsumption)
      return;

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
    if (!kPruneUnit)
      break;

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
    if (active[i] && (!kPruneUsefulness || useful[i]))
      E.push_back(inputE[i]);

  const ui eSize = (ui)E.size();
  const ui hSize = (ui)hitSets.size();
  addCliqueCountOrThrow(nontrivialSeedSolverCalls, 1);
  maximumSeedConstraints = max(maximumSeedConstraints, hSize);
#if !defined(PURE_HITSET_DYNAMIC)
  if (hSize > kHitsetCapacity) {
    if (usePivotFallback != nullptr) {
      *usePivotFallback = true;
      addCliqueCountOrThrow(solverCapacityFallbacks, 1);
      return {};
    }
    vector<vector<ui>> residual = backtrackingBranchBound(E, hitSets);
    for (vector<ui> &solution : residual) {
      solution.insert(solution.end(), preForced.begin(), preForced.end());
      sort(solution.begin(), solution.end());
    }
    return residual;
  }
#endif
#if !defined(PURE_HITSET_DYNAMIC)
  static constexpr ui fixedMaskWords = kHitsetCapacity / 64;
  using Mask = array<ull, fixedMaskWords>;
  const ui maskWords = fixedMaskWords;
  auto zeroMask = []() -> Mask { return Mask{}; };
#else
  using Mask = vector<ull>;
  const ui maskWords = (hSize + 63) >> 6;
  auto zeroMask = [&]() -> Mask { return Mask(maskWords, 0ULL); };
#endif

  auto setBit = [](Mask &m, ui bit) {
    m[bit >> 6] |= (1ULL << (bit & 63));
  };
  auto clearBit = [](Mask &m, ui bit) {
    m[bit >> 6] &= ~(1ULL << (bit & 63));
  };
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

  // cov[i] = bitmask of hitSets that E[i] covers.
  // Non-const so subsumption can clear implied constraint bits.
  Mask fullMask = zeroMask();
  for (ui bit = 0; bit < hSize; bit++)
    setBit(fullMask, bit);
  vector<Mask> cov(eSize, zeroMask());
  {
    if (++eIndexToken == 0) {
      fill(eIndexStamp.begin(), eIndexStamp.end(), 0);
      eIndexToken = 1;
    }
    for (ui i = 0; i < eSize; i++) {
      eIndex[E[i]] = i;
      eIndexStamp[E[i]] = eIndexToken;
    }
    for (ui bit = 0; bit < hSize; bit++) {
      // normalizeAndSubsume already leaves hitSets in size-then-lexicographic
      // order, so lower bits retain the same fail-first priority directly.
      for (ui v : hitSets[bit]) {
        if (eIndexStamp[v] == eIndexToken)
          setBit(cov[eIndex[v]], bit);
      }
    }
  }

  // Initial candidate order: descending coverage count so high-utility
  // vertices are tried first, finding solutions sooner for better pruning.
  vector<ui> initCands(eSize);
  iota(initCands.begin(), initCands.end(), 0);
  sort(initCands.begin(), initCands.end(), [&](ui a, ui b) {
    const int aCoverage = popcountMask(cov[a]);
    const int bCoverage = popcountMask(cov[b]);
    if (aCoverage != bCoverage)
      return aCoverage > bCoverage;
    return E[a] < E[b];
  });

  // ── Unit propagation ──────────────────────────────────────────────────────
  // For any constraint covered by exactly one candidate, that candidate is
  // forced into every solution. Pre-select all forced candidates, filter the
  // remaining candidates to be clique-compatible with them, and start the DFS
  // with the forced coverage already accumulated.
  vector<ui> forcedIdxs; // E-indices forced into every solution
  Mask forcedCov = zeroMask();

  if (kPruneUnit) {
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
            if (!adj(E[fv], E[sole])) {
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
            if (ci != sole && adj(E[sole], E[ci]))
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

  // ── Constraint subsumption ────────────────────────────────────────────────
  // Constraint i is subsumed by constraint j when every candidate covering j
  // also covers i (coverSet(j) ⊆ coverSet(i)). Satisfying j then implies
  // satisfying i, so i can be dropped from fullMask.
  // Work on the post-unit-propagation candidates so forced coverage is visible.
  if (kPruneSubsumption) {
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
        // subsumption
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

  // If forced coverage already satisfies all remaining constraints, return it.
  if (coversAll(forcedCov, fullMask)) {
    if (forcedIdxs.empty() && preForced.empty())
      return {}; // nothing forced, nothing to cover
    vector<ui> sol = preForced;
    for (ui idx : forcedIdxs)
      sol.push_back(E[idx]);
    sort(sol.begin(), sol.end());
    return {sol};
  }

  // Avoid building the quadratic compatibility table when the budgeted DFS
  // is already guaranteed to fall back at its root.
  if (solverWorkBudgetEnabled && usePivotFallback != nullptr) {
    const ull candidateCount = static_cast<ull>(initCands.size());
    if (candidateCount > 1 &&
        candidateCount * (candidateCount - 1) / 2 > solverWorkBudget) {
      const Mask uncovered = andNot(fullMask, forcedCov);
      Mask available = forcedCov;
      ull rootCandidates = 0;
      for (ui candidate : initCands) {
        orEq(available, cov[candidate]);
        if (!kPruneZeroCoverage || intersects(cov[candidate], uncovered))
          ++rootCandidates;
      }
      if (!kPruneFailFirst || coversAll(available, fullMask)) {
        const ull unavoidableRootWork =
            rootCandidates * (rootCandidates - 1) / 2;
        if (unavoidableRootWork > solverWorkBudget) {
          *usePivotFallback = true;
          if constexpr (kDiagnosticsEnabled) {
            addCliqueCountOrThrow(solverBudgetFallbacks, 1);
            addCliqueCountOrThrow(solverCertifiedBudgetFallbacks, 1);
          }
          return {};
        }
      } else {
        return {};
      }
    }
  }

  // compat[i*eSize+j] = 1 iff E[i] and E[j] are adjacent in the graph.
  vector<char> compat(static_cast<size_t>(eSize) * eSize, 0);
  for (ui i = 0; i < eSize; ++i) {
    const AdjacencyRow row = adjacentVertices(E[i]);
    for (ui j = i + 1; j < eSize; ++j)
      if (binary_search(row.begin(), row.end(), E[j]))
        compat[static_cast<size_t>(i) * eSize + j] =
            compat[static_cast<size_t>(j) * eSize + i] = 1;
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
      // Before recording, verify cur is not a superset of an existing solution.
      if (kPruneAntichain) {
        for (const auto &s : solutions)
          if (includes(cur.begin(), cur.end(), s.begin(), s.end()))
            return;
      }
      // Remove any existing solutions that cur dominates (cur is a subset).
      if (kPruneAntichain)
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

    // Superset pruning: cur already contains a known minimal solution so any
    // extension of cur cannot be minimal.
    if (kPruneAntichain) {
      for (const auto &s : solutions)
        if (includes(cur.begin(), cur.end(), s.begin(), s.end()))
          return;
    }

    // Fail-first dead-branch check: for every uncovered constraint verify at
    // least one candidate can cover it. If any constraint is impossible, prune.
    if (kPruneFailFirst) {
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
      if (kPruneZeroCoverage && !intersects(cov[ci], uncovered))
        continue;

      // Build next-level candidates: those in cands with E-index > ci that
      // are adjacent to ci (enforces clique property and avoids duplicates).
      vector<ui> &next = candScratch[depth];
      next.clear();
      if (next.capacity() < cands.size())
        next.reserve(cands.size());
      for (ui cj : cands)
        if (cj > ci) {
          if (solverWorkBudgetEnabled && usePivotFallback != nullptr &&
              solverWork >= solverWorkBudget) {
            budgetExceeded = true;
            break;
          }
          ++solverWork;
          if (compat[ci * eSize + cj])
            next.push_back(cj);
        }
      if (budgetExceeded)
        return;
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
    if constexpr (kDiagnosticsEnabled)
      addCliqueCountOrThrow(solverBudgetFallbacks, 1);
    return {};
  }

  if (!kPruneAntichain)
    solutions = minimalByInclusion(std::move(solutions));
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
  const AdjacencyRow row = adjacentVertices(u);
  // A short contiguous CSR scan is cheaper than one hash/binary lookup per P
  // vertex even when the row is moderately larger than P.
  if (row.size() <= P.size() * 32) {
    for (ui v : row)
      score += eIndexStamp[v] == eIndexToken;
  } else {
    for (ui v : P)
      score += adj(u, v);
  }
  return score;
}

void ReorderSib::scanPurePXRState(
    const vector<ui> &P, const vector<ui> &X, ui &pivot, ui &minPScore,
    ui &universalP, bool &xUniversal) {
  const ui pSize = static_cast<ui>(P.size());
  if (++eIndexToken == 0) {
    fill(eIndexStamp.begin(), eIndexStamp.end(), 0);
    eIndexToken = 1;
  }
  for (ui v : P)
    eIndexStamp[v] = eIndexToken;

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

void ReorderSib::pureMatchingParts(
    const vector<ui> &P, vector<ui> &forced,
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

// Exact constant-size kernel for a formal branch B=(M,Q), |Q| <= 4.  It
// enumerates every clique-compatible subset of Q and accepts it only when no
// graph vertex extends M union S.  This is the same global maximality
// condition enforced by the X set in ordinary PXR, without constructing P/X
// child vectors and recursive frames for at most sixteen possibilities.
void ReorderSib::enumerateSmallPureBranch(const vector<ui> &M,
                                          const vector<ui> &Q) {
  if (M.empty() || Q.size() > 4)
    throw logic_error("small Pure branch kernel received invalid dimensions");

  const AdjacencyRow firstRow = adjacentVertices(M[0]);
  vector<ui> common(firstRow.begin(), firstRow.end());
  vector<ui> scratch;
  for (size_t i = 1; i < M.size() && !common.empty(); ++i) {
    intersectInto(scratch, common, adjacentVertices(M[i]));
    common.swap(scratch);
  }

  const unsigned subsetCount = 1U << static_cast<unsigned>(Q.size());
  vector<ui> selected;
  vector<ui> clique;
  selected.reserve(Q.size());
  clique.reserve(M.size() + Q.size());
  for (unsigned mask = 0; mask < subsetCount; ++mask) {
    const size_t selectedCount =
        static_cast<size_t>(__builtin_popcount(mask));
    if (M.size() + selectedCount < minCliqueSize)
      continue;

    selected.clear();
    bool isClique = true;
    for (size_t i = 0; i < Q.size() && isClique; ++i) {
      if ((mask & (1U << i)) == 0)
        continue;
      for (ui chosen : selected) {
        if (!adj(chosen, Q[i])) {
          isClique = false;
          break;
        }
      }
      selected.push_back(Q[i]);
    }
    if (!isClique)
      continue;

    bool hasExtension = false;
    for (ui candidate : common) {
      bool extends = true;
      for (ui chosen : selected) {
        if (!adj(candidate, chosen)) {
          extends = false;
          break;
        }
      }
      if (extends) {
        hasExtension = true;
        break;
      }
    }
    if (hasExtension)
      continue;

    clique.assign(M.begin(), M.end());
    clique.insert(clique.end(), selected.begin(), selected.end());
    recordPureClique(clique);
  }
}

// Stop-after-one Bron--Kerbosch for a formal branch B=(M,Q). The initial X is
// the set of common neighbors of M that lie outside Q. Carrying X through the
// recursion makes every returned leaf globally maximal, not merely maximal
// inside M union Q.
bool ReorderSib::findOnePure(const vector<ui> &M, const vector<ui> &Q,
                             vector<ui> &found) {
  found.clear();
  if (M.empty())
    return false;

  const AdjacencyRow firstRow = adjacentVertices(M[0]);
  vector<ui> common(firstRow.begin(), firstRow.end());
  vector<ui> scratch;
  for (ui i = 1; i < (ui)M.size() && !common.empty(); i++) {
    intersectInto(scratch, common, adjacentVertices(M[i]));
    common.swap(scratch);
  }

  vector<ui> X;
  setDiffInto(X, common, Q);
  vector<ui> R = M;
  vector<ui> P = Q;
  const size_t bufferCount = P.size() + 1;
  if (pxrPBuffers.size() < bufferCount) {
    pxrPBuffers.resize(bufferCount);
    pxrXBuffers.resize(bufferCount);
  }
  return findOnePureRecursive(R, P, X, found, 0);
}

bool ReorderSib::findOnePureRecursive(vector<ui> &R, vector<ui> &P,
                                      vector<ui> &X, vector<ui> &found,
                                      size_t depth) {
  incrementSearchStateOrThrow(checksCount);
  incrementSearchStateOrThrow(findOnePxrStates);
  if (R.size() + P.size() < minCliqueSize)
    return false;

  // These structural terminals are enabled by default in the Pure PXR lane
  // and do not depend on the Hybrid graph-level portfolio. The master ET
  // ablation used by the controlled PXR-state experiment disables them.
  if (P.empty()) {
    if (X.empty()) {
      found = R;
      sort(found.begin(), found.end());
      return true;
    }
    return false;
  }

  // A zero/one-candidate child can be decided without another BK level.
  if (kEt1Enabled && P.size() == 1) {
    const ui extension = P.front();
    for (ui x : X)
      if (adj(x, extension))
        return false;
    incrementSearchStateOrThrow(et1EnumeratedStates);
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
  if (kEt1Enabled && minPScore + 1 == pSize) {
    incrementSearchStateOrThrow(et1EnumeratedStates);
    found = R;
    found.insert(found.end(), P.begin(), P.end());
    sort(found.begin(), found.end());
    return found.size() >= minCliqueSize;
  }

  // The complement of P is a matching: isolated complement vertices are
  // forced and one endpoint from every missing edge gives a maximal clique.
  if (kEt2Enabled && X.empty() && minPScore + 2 >= pSize) {
    if constexpr (kDiagnosticsEnabled)
      ++findOne2PlexTerminals;
    vector<ui> forced;
    vector<pair<ui, ui>> missingEdges;
    pureMatchingParts(P, forced, missingEdges);
    incrementSearchStateOrThrow(et2EnumeratedStates);
    found = R;
    found.insert(found.end(), forced.begin(), forced.end());
    for (const auto &edge : missingEdges)
      found.push_back(edge.first);
    sort(found.begin(), found.end());
    return found.size() >= minCliqueSize;
  }

  // If every P vertex misses at most two P-neighbors, the complement consists
  // of paths and cycles. Reuse the exact 3-plex DP to obtain one witness.
  if (kEt3Enabled && X.empty() && minPScore + 3 >= pSize &&
      (!kEt2Enabled || minPScore + 2 < pSize)) {
    FastPlex3Result plex = solveFastPlex3Subtree(
        adjVertices, adjOffsets, adjHash, P, static_cast<ui>(R.size()), &R,
        minCliqueSize, nullptr);
    if (plex.handled) {
      if constexpr (kDiagnosticsEnabled)
        ++findOne3PlexTerminals;
      found = std::move(plex.witness);
      if (plex.found)
        incrementSearchStateOrThrow(et3EnumeratedStates);
      sort(found.begin(), found.end());
      return plex.found;
    }
  }

  // A P-universal vertex must be present in every maximal continuation, so
  // force it into R and make a single recursive call.
  if (universalP != numeric_limits<ui>::max()) {
    vector<ui> &childP = pxrPBuffers[depth + 1];
    vector<ui> &childX = pxrXBuffers[depth + 1];
    intersectInto(childP, P, adjacentVertices(universalP));
    intersectInto(childX, X, adjacentVertices(universalP));
    R.push_back(universalP);
    const bool result =
        findOnePureRecursive(R, childP, childX, found, depth + 1);
    R.pop_back();
    return result;
  }

  vector<ui> branchRoots;
  setDiffInto(branchRoots, P, adjacentVertices(pivot));
  for (ui v : branchRoots) {
    vector<ui> &childP = pxrPBuffers[depth + 1];
    vector<ui> &childX = pxrXBuffers[depth + 1];
    intersectInto(childP, P, adjacentVertices(v));
    intersectInto(childX, X, adjacentVertices(v));
    R.push_back(v);
    if (findOnePureRecursive(R, childP, childX, found, depth + 1)) {
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

  const AdjacencyRow firstRow = adjacentVertices(M[0]);
  vector<ui> common(firstRow.begin(), firstRow.end());
  vector<ui> scratch;
  for (ui i = 1; i < (ui)M.size() && !common.empty(); i++) {
    intersectInto(scratch, common, adjacentVertices(M[i]));
    common.swap(scratch);
  }

  vector<ui> X;
  setDiffInto(X, common, Q);
  vector<ui> R = M;
  vector<ui> P = Q;
  const size_t bufferCount = P.size() + 1;
  if (pxrPBuffers.size() < bufferCount) {
    pxrPBuffers.resize(bufferCount);
    pxrXBuffers.resize(bufferCount);
  }
  enumerateAllPureBranchRecursive(R, P, X, 0);
}

void ReorderSib::enumerateAllPureBranchRecursive(
    vector<ui> &R, vector<ui> &P, vector<ui> &X, size_t depth) {
  incrementSearchStateOrThrow(checksCount);
  incrementSearchStateOrThrow(fullPxrStates);
  if (R.size() + P.size() < minCliqueSize)
    return;

  // The same default terminals used by witness search also consume an
  // over-capacity Pure branch without descending through ordinary Pivot-BK.
  if (P.empty()) {
    if (X.empty())
      recordPureClique(R);
    return;
  }

  if (kEt1Enabled && P.size() == 1) {
    const ui extension = P.front();
    for (ui x : X)
      if (adj(x, extension))
        return;
    incrementSearchStateOrThrow(et1EnumeratedStates);
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

  if (kEt1Enabled && minPScore + 1 == pSize) {
    incrementSearchStateOrThrow(et1EnumeratedStates);
    vector<ui> clique = R;
    clique.insert(clique.end(), P.begin(), P.end());
    if (clique.size() >= minCliqueSize)
      recordPureClique(std::move(clique));
    return;
  }

  if (kEt2Enabled && X.empty() && minPScore + 2 >= pSize) {
    if constexpr (kDiagnosticsEnabled)
      ++fullPxr2PlexTerminals;
    vector<ui> forced;
    vector<pair<ui, ui>> missingEdges;
    pureMatchingParts(P, forced, missingEdges);
    vector<ui> clique = R;
    clique.insert(clique.end(), forced.begin(), forced.end());
    function<void(size_t)> materialize = [&](size_t at) {
      if (at == missingEdges.size()) {
        incrementSearchStateOrThrow(et2EnumeratedStates);
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

  if (kEt3Enabled && X.empty() && minPScore + 3 >= pSize &&
      (!kEt2Enabled || minPScore + 2 < pSize)) {
    FastCliqueSink sink =
        [&](const vector<ui> &clique) { recordPureClique(clique); };
    FastPlex3Result plex = solveFastPlex3Subtree(
        adjVertices, adjOffsets, adjHash, P, static_cast<ui>(R.size()), &R,
        minCliqueSize, &sink);
    if (plex.handled) {
      if constexpr (kDiagnosticsEnabled)
        ++fullPxr3PlexTerminals;
      addSearchStatesOrThrow(et3EnumeratedStates,
                             plex.enumeratedCliqueCount);
      return;
    }
  }

  if (universalP != numeric_limits<ui>::max()) {
    vector<ui> &childP = pxrPBuffers[depth + 1];
    vector<ui> &childX = pxrXBuffers[depth + 1];
    intersectInto(childP, P, adjacentVertices(universalP));
    intersectInto(childX, X, adjacentVertices(universalP));
    R.push_back(universalP);
    enumerateAllPureBranchRecursive(R, childP, childX, depth + 1);
    R.pop_back();
    return;
  }

  vector<ui> branchRoots;
  setDiffInto(branchRoots, P, adjacentVertices(pivot));
  for (ui v : branchRoots) {
    vector<ui> &childP = pxrPBuffers[depth + 1];
    vector<ui> &childX = pxrXBuffers[depth + 1];
    intersectInto(childP, P, adjacentVertices(v));
    intersectInto(childX, X, adjacentVertices(v));
    R.push_back(v);
    enumerateAllPureBranchRecursive(R, childP, childX, depth + 1);
    R.pop_back();

    auto pIt = lower_bound(P.begin(), P.end(), v);
    if (pIt != P.end() && *pIt == v)
      P.erase(pIt);
    X.insert(lower_bound(X.begin(), X.end(), v), v);
  }
}

static ull hashClique(const vector<ui> &clique) {
  ull hash = 1469598103934665603ULL;
  for (ui vertex : clique) {
    hash ^= static_cast<ull>(vertex);
    hash *= 1099511628211ULL;
  }
  hash ^= static_cast<ull>(clique.size());
  hash *= 1099511628211ULL;
  return hash;
}

static ull mixCliqueHash(ull hash) {
  hash ^= hash >> 30;
  hash *= 0xbf58476d1ce4e5b9ULL;
  hash ^= hash >> 27;
  hash *= 0x94d049bb133111ebULL;
  hash ^= hash >> 31;
  return hash;
}

size_t ReorderSib::emittedCliqueHashSlot(ull hash) const {
  const size_t mask = emittedHashHeads.size() - 1;
  size_t slot = static_cast<size_t>(mixCliqueHash(hash)) & mask;
  while (emittedHashHeads[slot] != numeric_limits<size_t>::max() &&
         emittedHashKeys[slot] != hash)
    slot = (slot + 1) & mask;
  return slot;
}

void ReorderSib::rehashEmittedCliqueIndex(size_t capacity) {
  vector<ull> oldKeys = std::move(emittedHashKeys);
  vector<size_t> oldHeads = std::move(emittedHashHeads);
  emittedHashKeys.assign(capacity, 0);
  emittedHashHeads.assign(capacity, numeric_limits<size_t>::max());
  emittedHashSlotsUsed = 0;

  for (size_t oldSlot = 0; oldSlot < oldHeads.size(); ++oldSlot) {
    if (oldHeads[oldSlot] == numeric_limits<size_t>::max())
      continue;
    const size_t slot = emittedCliqueHashSlot(oldKeys[oldSlot]);
    emittedHashKeys[slot] = oldKeys[oldSlot];
    emittedHashHeads[slot] = oldHeads[oldSlot];
    ++emittedHashSlotsUsed;
  }
}

bool ReorderSib::recordPureClique(vector<ui> C) {
  sort(C.begin(), C.end());
  const ull hash = hashClique(C);

  if (emittedHashHeads.empty())
    rehashEmittedCliqueIndex(16);
  size_t slot = emittedCliqueHashSlot(hash);
  if (emittedHashHeads[slot] != numeric_limits<size_t>::max()) {
    for (size_t cliqueId = emittedHashHeads[slot];
         cliqueId != numeric_limits<size_t>::max();
         cliqueId = emittedHashNext[cliqueId]) {
      if (!storedCliqueEquals(static_cast<ui>(cliqueId), C))
        continue;
      if constexpr (kDiagnosticsEnabled)
        addCliqueCountOrThrow(dupBlocked, 1);
      return false;
    }
  }

  if (storedCliqueCount() > numeric_limits<ui>::max())
    throw overflow_error("materialized clique index exceeds uint32_t");

  const ui cliqueId = static_cast<ui>(storedCliqueCount());
  addCliqueCountOrThrow(cliqueCount, 1);
  if constexpr (kDiagnosticsEnabled)
    maxCliqueSize = max(maxCliqueSize, C.size());
  cliqueVertices.insert(cliqueVertices.end(), C.begin(), C.end());
  cliqueOffsets.push_back(cliqueVertices.size());

  if (emittedHashHeads[slot] == numeric_limits<size_t>::max() &&
      (emittedHashSlotsUsed + 1) * 10 > emittedHashHeads.size() * 7) {
    rehashEmittedCliqueIndex(emittedHashHeads.size() * 2);
    slot = emittedCliqueHashSlot(hash);
  }
  const size_t previousHead = emittedHashHeads[slot];
  emittedHashNext.push_back(previousHead);
  if (previousHead == numeric_limits<size_t>::max()) {
    emittedHashKeys[slot] = hash;
    ++emittedHashSlotsUsed;
  }
  emittedHashHeads[slot] = cliqueId;
  for (ui v : C)
    cliqueIdsByVertex[v].push_back(cliqueId);
  return true;
}

vector<vector<ui>> ReorderSib::getCliques() const {
  vector<vector<ui>> restored;
  restored.reserve(storedCliqueCount());
  for (size_t cliqueId = 0; cliqueId < storedCliqueCount(); ++cliqueId) {
    const auto begin = cliqueVertices.begin() + cliqueOffsets[cliqueId];
    const auto end = cliqueVertices.begin() + cliqueOffsets[cliqueId + 1];
    vector<ui> clique(begin, end);
    for (ui &vertex : clique)
      vertex = internalToOriginal[vertex];
    sort(clique.begin(), clique.end());
    restored.push_back(std::move(clique));
  }
  return restored;
}

void ReorderSib::findAllMaximalCliquesPure() {
  cliqueCount = 0;
  dupBlocked = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  findOnePxrStates = 0;
  fullPxrStates = 0;
  et1EnumeratedStates = 0;
  et2EnumeratedStates = 0;
  et3EnumeratedStates = 0;
  solverBudgetFallbacks = 0;
  solverCertifiedBudgetFallbacks = 0;
  solverCapacityFallbacks = 0;
  nontrivialSeedSolverCalls = 0;
  maximumSeedConstraints = 0;
  findOne2PlexTerminals = 0;
  findOne3PlexTerminals = 0;
  fullPxr2PlexTerminals = 0;
  fullPxr3PlexTerminals = 0;
  worklistPushes = 0;
  worklistPops = 0;
  maximumWorklistSize = 0;
  minSizePrunedBranches = 0;
  coverLookupCalls = 0;
  findOneCalls = 0;
  findOneSuccesses = 0;
  findOneCliqueSizeTotal = 0;
  maximumFindOneCliqueSize = 0;
  seedSolverCalls = 0;
  generatedBranches = 0;
  fullPxrFallbackBranches = 0;
  smallQFullPxrBranches = 0;
  poppedQSizeBuckets.fill(0);
  cliqueVertices.clear();
  cliqueOffsets.clear();
  cliqueOffsets.push_back(0);
  emittedHashKeys.clear();
  emittedHashHeads.clear();
  emittedHashNext.clear();
  emittedHashSlotsUsed = 0;
  cliqueIdsByVertex.assign(n, {});

  vector<PureBranch> worklist;
  vector<PureBranch> freeBranches;
  worklist.reserve(256);
  freeBranches.reserve(256);
  auto acquireBranch = [&]() {
    if (freeBranches.empty())
      return PureBranch{};
    PureBranch branch = std::move(freeBranches.back());
    freeBranches.pop_back();
    return branch;
  };
  auto recycleBranch = [&](PureBranch branch) {
    branch.mustin.clear();
    branch.expandTo.clear();
    freeBranches.push_back(std::move(branch));
  };
  auto pushBranch = [&](PureBranch nextBranch) {
    incrementSearchStateOrThrow(worklistPushes);
    worklist.push_back(std::move(nextBranch));
    if constexpr (kDiagnosticsEnabled)
      maximumWorklistSize = max(maximumWorklistSize, worklist.size());
  };

#if !defined(PURE_LEAN_BENCHMARK)
  auto t0 = chrono::high_resolution_clock::now();
#endif
  // The old implementation materialized all n root branches up front. A LIFO
  // stack always completed root v's descendants before visiting root v + 1,
  // so creating one root at a time preserves the exact traversal while
  // avoiding O(n) live vector objects and one copied forward-neighbor list per
  // unvisited root.
  for (ui root = 0; root < n; ++root) {
    const AdjacencyRow neighbors = adjacentVertices(root);
    const auto forwardBegin = neighbors.begin() + firstForwardNeighbor[root];
    const size_t forwardCount =
        static_cast<size_t>(neighbors.end() - forwardBegin);
    if (1 + forwardCount < minCliqueSize) {
      incrementSearchStateOrThrow(minSizePrunedBranches);
      continue;
    }
    PureBranch rootBranch = acquireBranch();
    rootBranch.mustin.clear();
    rootBranch.expandTo.clear();
    rootBranch.mustin.push_back(root);
    rootBranch.expandTo.assign(forwardBegin, neighbors.end());
    pushBranch(std::move(rootBranch));
    while (!worklist.empty()) {
      PureBranch branch = std::move(worklist.back());
      worklist.pop_back();
      incrementSearchStateOrThrow(worklistPops);
      const size_t qBucket = min<size_t>(branch.expandTo.size(), 7);
      incrementSearchStateOrThrow(poppedQSizeBuckets[qBucket]);

      if (branch.mustin.size() + branch.expandTo.size() < minCliqueSize) {
        incrementSearchStateOrThrow(minSizePrunedBranches);
        recycleBranch(std::move(branch));
        continue;
      }

      if (smallQFullPxrThreshold != 0 &&
          branch.expandTo.size() <= smallQFullPxrThreshold) {
        incrementSearchStateOrThrow(smallQFullPxrBranches);
        if (branch.expandTo.size() <= 4)
          enumerateSmallPureBranch(branch.mustin, branch.expandTo);
        else
          enumerateAllPureBranch(branch.mustin, branch.expandTo);
        recycleBranch(std::move(branch));
        continue;
      }

      incrementSearchStateOrThrow(coverLookupCalls);
      vector<ui> covers = collectAllCoveringCliques(branch.mustin);
      if (covers.empty()) {
        vector<ui> found;
        incrementSearchStateOrThrow(findOneCalls);
        if (findOnePure(branch.mustin, branch.expandTo, found)) {
          incrementSearchStateOrThrow(findOneSuccesses);
          addSearchStatesOrThrow(findOneCliqueSizeTotal,
                                 static_cast<ull>(found.size()));
          if constexpr (kDiagnosticsEnabled)
            maximumFindOneCliqueSize =
                max(maximumFindOneCliqueSize, found.size());
          recordPureClique(std::move(found));
          // The unchanged branch retains every unseen target. On its next pop,
          // the clique just recorded necessarily covers its must-in set.
          pushBranch(std::move(branch));
        } else {
          recycleBranch(std::move(branch));
        }
        continue;
      }

      bool usePivotFallback = false;
      incrementSearchStateOrThrow(seedSolverCalls);
      vector<vector<ui>> seeds = generateExactSiblingSets(
          branch.expandTo, covers, &usePivotFallback);
      if (usePivotFallback) {
        incrementSearchStateOrThrow(fullPxrFallbackBranches);
        enumerateAllPureBranch(branch.mustin, branch.expandTo);
        recycleBranch(std::move(branch));
        continue;
      }
      // Reverse insertion preserves the solver's deterministic seed order
      // under the LIFO worklist; correctness does not depend on this order.
      for (auto it = seeds.rbegin(); it != seeds.rend(); ++it) {
        PureBranch nextBranch = acquireBranch();
        unionSetInto(nextBranch.mustin, branch.mustin, *it);
        commonExpandInto(nextBranch.expandTo, branch.expandTo, *it);
        incrementSearchStateOrThrow(generatedBranches);
        if (nextBranch.mustin.size() + nextBranch.expandTo.size() <
            minCliqueSize) {
          incrementSearchStateOrThrow(minSizePrunedBranches);
          recycleBranch(std::move(nextBranch));
          continue;
        }
        pushBranch(std::move(nextBranch));
      }
      recycleBranch(std::move(branch));
    }
  }
#if !defined(PURE_LEAN_BENCHMARK)
  auto t1 = chrono::high_resolution_clock::now();
  const double ms = chrono::duration<double, milli>(t1 - t0).count();
  // PXR states are recursive (R,P,X) entries. ET states are terminal clique
  // possibilities directly enumerated after an early-termination hit.
  ull pxrStates = 0;
  addSearchStatesOrThrow(pxrStates, findOnePxrStates);
  addSearchStatesOrThrow(pxrStates, fullPxrStates);
  if (pxrStates != checksCount)
    throw logic_error("PXR state counters disagree with checksCount");
  ull etStates = 0;
  addSearchStatesOrThrow(etStates, et1EnumeratedStates);
  addSearchStatesOrThrow(etStates, et2EnumeratedStates);
  addSearchStatesOrThrow(etStates, et3EnumeratedStates);
  ull totalStates = pxrStates;
  addSearchStatesOrThrow(totalStates, etStates);
#endif

  if (static_cast<ull>(storedCliqueCount()) != cliqueCount)
    throw logic_error("stored Pure clique count differs from numeric count");

#if defined(PURE_LEAN_BENCHMARK)
  cout << "pure.cliques=" << cliqueCount << '\n'
       << "pure.stored_cliques=" << storedCliqueCount() << '\n'
       << "pure.minimum_clique_size=" << minCliqueSize << '\n'
       << "pure.config.budget=";
  if (solverWorkBudgetEnabled)
    cout << solverWorkBudget;
  else
    cout << "unlimited";
  cout << '\n' << "pure.config.lean_benchmark=1" << endl;
#else
  cout << fixed << setprecision(3) << "PureReorderSib: cliques=" << cliqueCount
       << "  storedCliques=" << storedCliqueCount()
       << "  dups=" << dupBlocked << "  maxSize=" << maxCliqueSize
       << "  minSize=" << minCliqueSize << "  checks=" << checksCount
       << "  budgetFallbacks=" << solverBudgetFallbacks
       << "  capacityFallbacks=" << solverCapacityFallbacks
       << "  time=" << ms << " ms" << endl;

  cout << "pure.cliques=" << cliqueCount << '\n'
       << "pure.stored_cliques=" << storedCliqueCount() << '\n'
       << "pure.duplicates_blocked=" << dupBlocked << '\n'
       << "pure.maximum_clique_size=" << maxCliqueSize << '\n'
       << "pure.minimum_clique_size=" << minCliqueSize << '\n'
       << "pure.checks=" << checksCount << '\n'
       << "pure.budget_fallbacks=" << solverBudgetFallbacks << '\n'
       << "pure.certified_budget_fallbacks="
       << solverCertifiedBudgetFallbacks << '\n'
       << "pure.capacity_fallbacks=" << solverCapacityFallbacks << '\n'
       << "pure.seed_solver_calls=" << nontrivialSeedSolverCalls << '\n'
       << "pure.maximum_seed_constraints=" << maximumSeedConstraints << '\n'
       << "pure.runtime_ms=" << ms << '\n'
       << "pure.config.method="
       << (method == SibMethod::OPTIMIZED ? "optimized" : "backtracking")
       << '\n'
       << "pure.config.budget=";
  if (solverWorkBudgetEnabled)
    cout << solverWorkBudget;
  else
    cout << "unlimited";
  cout << '\n'
       << "pure.config.small_q_full_pxr_threshold="
       << smallQFullPxrThreshold << '\n'
       << "pure.config.hitset_capacity=";
  if (kHitsetDynamic)
    cout << "dynamic";
  else
    cout << kHitsetCapacity;
  cout << '\n'
       << "pure.config.hitset_dynamic=" << kHitsetDynamic << '\n'
       << "pure.config.et=" << (kEt1Enabled || kEt2Enabled || kEt3Enabled)
       << '\n'
       << "pure.config.et1=" << kEt1Enabled << '\n'
       << "pure.config.et2=" << kEt2Enabled << '\n'
       << "pure.config.et3=" << kEt3Enabled << '\n'
       << "pure.config.et.findone_2plex=" << kEt2Enabled << '\n'
       << "pure.config.et.findone_3plex=" << kEt3Enabled << '\n'
       << "pure.config.et.full_pxr_2plex=" << kEt2Enabled << '\n'
       << "pure.config.et.full_pxr_3plex=" << kEt3Enabled << '\n'
       << "pure.config.pruning.normalization=" << kPruneNormalization << '\n'
       << "pure.config.pruning.subsumption=" << kPruneSubsumption << '\n'
       << "pure.config.pruning.unit=" << kPruneUnit << '\n'
       << "pure.config.pruning.usefulness=" << kPruneUsefulness << '\n'
       << "pure.config.pruning.antichain=" << kPruneAntichain << '\n'
       << "pure.config.pruning.fail_first=" << kPruneFailFirst << '\n'
       << "pure.config.pruning.zero_coverage=" << kPruneZeroCoverage << '\n'
       << "pure.counter.findone_2plex_terminals=" << findOne2PlexTerminals
       << '\n'
       << "pure.counter.findone_3plex_terminals=" << findOne3PlexTerminals
       << '\n'
       << "pure.counter.full_pxr_2plex_terminals=" << fullPxr2PlexTerminals
       << '\n'
       << "pure.counter.full_pxr_3plex_terminals=" << fullPxr3PlexTerminals
       << '\n'
       << "pure.counter.worklist_pushes=" << worklistPushes << '\n'
       << "pure.counter.worklist_pops=" << worklistPops << '\n'
       << "pure.counter.worklist_max_size=" << maximumWorklistSize << '\n'
       << "pure.counter.worklist_min_size_pruned=" << minSizePrunedBranches << '\n'
       << "pure.counter.cover_lookups=" << coverLookupCalls << '\n'
       << "pure.counter.findone_calls=" << findOneCalls << '\n'
       << "pure.counter.findone_successes=" << findOneSuccesses << '\n'
       << "pure.counter.findone_clique_size_total=" << findOneCliqueSizeTotal << '\n'
       << "pure.counter.findone_clique_size_max=" << maximumFindOneCliqueSize << '\n'
       << "pure.counter.seed_solver_calls=" << seedSolverCalls << '\n'
       << "pure.counter.generated_branches=" << generatedBranches << '\n'
       << "pure.counter.full_pxr_fallback_branches=" << fullPxrFallbackBranches << '\n'
       << "pure.counter.small_q_full_pxr_branches=" << smallQFullPxrBranches << '\n'
       << "pure.counter.worklist_q_size_0=" << poppedQSizeBuckets[0] << '\n'
       << "pure.counter.worklist_q_size_1=" << poppedQSizeBuckets[1] << '\n'
       << "pure.counter.worklist_q_size_2=" << poppedQSizeBuckets[2] << '\n'
       << "pure.counter.worklist_q_size_3=" << poppedQSizeBuckets[3] << '\n'
       << "pure.counter.worklist_q_size_4=" << poppedQSizeBuckets[4] << '\n'
       << "pure.counter.worklist_q_size_5=" << poppedQSizeBuckets[5] << '\n'
       << "pure.counter.worklist_q_size_6=" << poppedQSizeBuckets[6] << '\n'
       << "pure.counter.worklist_q_size_ge_7=" << poppedQSizeBuckets[7] << '\n'
       << "pure.state.findone_pxr=" << findOnePxrStates << '\n'
       << "pure.state.full_pxr=" << fullPxrStates << '\n'
       << "pure.state.pxr=" << pxrStates << '\n'
       << "pure.state.et1=" << et1EnumeratedStates << '\n'
       << "pure.state.et2=" << et2EnumeratedStates << '\n'
       << "pure.state.et3=" << et3EnumeratedStates << '\n'
       << "pure.state.et=" << etStates << '\n'
       << "pure.state.total=" << totalStates << endl;
#endif

}
