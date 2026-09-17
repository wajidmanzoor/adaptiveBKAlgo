#include "../inc/helpers.h"
#include "../inc/config.h"
#include "../inc/fast_plex3.h"
#include <functional>
#include <numeric>
#include <type_traits>

namespace {
using pure_config::kEt1Enabled;
using pure_config::kEt2Enabled;
using pure_config::kEt3Enabled;
using pure_config::kAdjHashThreshold;
using pure_config::kAdaptiveDirectMinCliquesPerRootDenominator;
using pure_config::kAdaptiveDirectQThreshold;
using pure_config::kAdaptiveDirectWarmupRoots;
using pure_config::kCoverCollectionCutoff;
using pure_config::kHitsetCapacity;
using pure_config::kHitsetDynamic;
using pure_config::kSmallQCcrThreshold;
using pure_config::kPruneAntichain;
using pure_config::kPruneFailFirst;
using pure_config::kPruneNormalization;
using pure_config::kPruneSubsumption;
using pure_config::kPruneUnit;
using pure_config::kPruneUsefulness;
using pure_config::kPruneZeroCoverage;

constexpr size_t kBinaryIntersectionRatio = 16;

void addCcrMetric(ull &metric, ull increment) {
  ull updated = 0;
  if (!tryAddUll(metric, increment, updated))
    throw overflow_error("CCRMCE experiment metric exceeds uint64_t");
  metric = updated;
}

} // namespace

bool ReorderSib::FlatAdjacencyHash::contains(ui vertex) const {
  constexpr ull kMultiplier = 11400714819323198485ULL;
  size_t slot = static_cast<size_t>(static_cast<ull>(vertex) * kMultiplier) &
                mask;
  while (slots[slot] != numeric_limits<ui>::max()) {
    if (slots[slot] == vertex)
      return true;
    slot = (slot + 1) & mask;
  }
  return false;
}

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
ReorderSib::ReorderSib(Graph &g, ui minCliqueSize)
    : minCliqueSize(max<ui>(1, minCliqueSize)) {
  n = g.n;
  cliqueCount = 0;
  solverWorkBudget = 0;
  solverWorkBudgetEnabled = false;
  solverBudgetFallbacks = 0;
  solverCapacityFallbacks = 0;
  seedSolverCalls = 0;
  maximumSeedConstraints = 0;
  ccrFindOneCalls = 0;
  ccrFullCalls = 0;
  ccrFindOneStates = 0;
  ccrFullStates = 0;
  ccrCoreExtractions = 0;
  ccrCoreVertices = 0;
  ccrResidualVertices = 0;

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
    if (row.size() > numeric_limits<size_t>::max() / 2)
      throw length_error("adjacency hash capacity overflow");
    const size_t needed = row.size() * 2;
    size_t capacity = 1;
    while (capacity < needed) {
      if (capacity > numeric_limits<size_t>::max() / 2)
        throw length_error("adjacency hash capacity overflow");
      capacity *= 2;
    }
    auto hash = make_unique<FlatAdjacencyHash>();
    hash->slots.assign(capacity, numeric_limits<ui>::max());
    hash->mask = capacity - 1;
    constexpr ull kMultiplier = 11400714819323198485ULL;
    for (ui vertex : row) {
      size_t slot =
          static_cast<size_t>(static_cast<ull>(vertex) * kMultiplier) &
          hash->mask;
      while (hash->slots[slot] != numeric_limits<ui>::max())
        slot = (slot + 1) & hash->mask;
      hash->slots[slot] = vertex;
    }
    adjHash[u] = std::move(hash);
  }

  eIndex.assign(n, 0);
  eIndexStamp.assign(n, 0);
  eIndexToken = 0;
  ccrIndex.assign(n, 0);
  ccrIndexStamp.assign(n, 0);
  ccrIndexToken = 0;
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
// weak-cover filters are applied here. Large exact constraint sets commonly
// end in the direct-CCRMCE fallback anyway. Stop once the configured cutoff is
// exceeded and let the caller take the same exact fallback immediately.
bool ReorderSib::collectAllCoveringCliques(const vector<ui> &M,
                                           vector<ui> &result) {
  result.clear();
  if (M.empty())
    return true;

  ui seed = n;
  for (ui v : M) {
    const size_t forwardCount =
        adjacentVertices(v).size() - firstForwardNeighbor[v];
    if (forwardCount <= kSmallQCcrThreshold)
      continue;
    if (seed == n ||
        cliqueIdsByVertex[v].size() < cliqueIdsByVertex[seed].size())
      seed = v;
  }
  if (seed == n)
    throw logic_error("cover lookup has no indexed branch vertex");

  const auto &posting = cliqueIdsByVertex[seed];
  const size_t neededCapacity = min(
      posting.size(), static_cast<size_t>(kCoverCollectionCutoff) + 1);
  if (result.capacity() < neededCapacity)
    result.reserve(neededCapacity);
  if (posting.size() <= kCoverCollectionCutoff) {
    for (ui cliqueId : posting)
      if (storedCliqueContains(cliqueId, M))
        result.push_back(cliqueId);
    return true;
  }
  for (ui cliqueId : posting) {
    if (storedCliqueContains(cliqueId, M)) {
      if (result.size() == kCoverCollectionCutoff)
        return false;
      result.push_back(cliqueId);
    }
  }
  return true;
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

  // The normal exact solver has a fixed 128-constraint capacity.  Build its
  // coverage masks directly when the unreduced problem already fits; larger
  // inputs retain the vector-set preprocessing because exact normalization
  // can still reduce them below the capacity.
  if (coveringCliqueIds.size() <= 64)
    return efficientHittingSetDirect(E, coveringCliqueIds,
                                     usePivotFallback);

  vector<vector<ui>> hitSets = buildHitSets(E, coveringCliqueIds);
  for (const vector<ui> &hitSet : hitSets)
    if (hitSet.empty())
      return {};

  if (coveringCliqueIds.size() == 1)
    return singletonBranches(hitSets[0]);
  return efficientHittingSet(E, std::move(hitSets), usePivotFallback);
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
  addCcrMetric(seedSolverCalls, 1);
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
          addCcrMetric(solverBudgetFallbacks, 1);
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
    addCcrMetric(solverBudgetFallbacks, 1);
    return {};
  }

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
// Optimizations used by the exact seed solver:
//   1. Bitmask coverage  — done-check and update are O(1) bitwise ops.
//   2. Incremental candidate list — compat-filtered frontier passed down,
//      no per-step binary_search into reordered adjacency.
//   3. Fail-first dead-branch  — prune as soon as any uncovered constraint
//      has zero candidates left.
//   4. "Covers nothing new" skip — a vertex that adds no new coverage can
//      never be part of a minimal solution; skip it unconditionally.
//   5. Live minimal-set maintenance — dominated solutions are removed as soon
//      as a smaller one is found.
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
  addCcrMetric(seedSolverCalls, 1);
  maximumSeedConstraints = max(maximumSeedConstraints, hSize);
#if !defined(PURE_HITSET_DYNAMIC)
  if (hSize > kHitsetCapacity) {
    if (usePivotFallback != nullptr) {
      *usePivotFallback = true;
      addCcrMetric(solverCapacityFallbacks, 1);
    }
    return {};
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
          addCcrMetric(solverBudgetFallbacks, 1);
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
    addCcrMetric(solverBudgetFallbacks, 1);
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

namespace {

ui ccrBitCount(const vector<ull> &bits) {
  ull count = 0;
  for (ull word : bits)
    count += static_cast<ull>(__builtin_popcountll(word));
  if (count > numeric_limits<ui>::max())
    throw overflow_error("CCRMCE local set exceeds uint32_t");
  return static_cast<ui>(count);
}

bool ccrBitIsSet(const vector<ull> &bits, ui vertex) {
  return (bits[vertex >> 6] & (1ULL << (vertex & 63))) != 0;
}

void ccrSetBit(vector<ull> &bits, ui vertex) {
  bits[vertex >> 6] |= 1ULL << (vertex & 63);
}

void ccrClearBit(vector<ull> &bits, ui vertex) {
  bits[vertex >> 6] &= ~(1ULL << (vertex & 63));
}

ui ccrFirstBit(const vector<ull> &bits) {
  for (size_t wordIndex = 0; wordIndex < bits.size(); ++wordIndex) {
    if (bits[wordIndex] == 0)
      continue;
    return static_cast<ui>(
        wordIndex * 64 + __builtin_ctzll(bits[wordIndex]));
  }
  throw logic_error("CCRMCE requested a vertex from an empty set");
}

template <typename Visitor>
void ccrForEachBit(const vector<ull> &bits, Visitor &&visit) {
  for (size_t wordIndex = 0; wordIndex < bits.size(); ++wordIndex) {
    ull word = bits[wordIndex];
    while (word != 0) {
      const unsigned bit = static_cast<unsigned>(__builtin_ctzll(word));
      visit(static_cast<ui>(wordIndex * 64 + bit));
      word &= word - 1;
    }
  }
}

} // namespace

// One-word CCRMCE for the high-frequency Q<=64 route. External X stays lazy:
// maximality normally rejects after very few adjacency tests, which is
// cheaper here than materializing every X-by-Q relation.
bool ReorderSib::runSmallCcrBranch(const vector<ui> &M,
                                   const vector<ui> &Q,
                                   const vector<ui> &common,
                                   vector<ui> *found) {
  if (Q.size() > kSmallQCcrThreshold)
    throw logic_error("one-word CCRMCE received too many candidates");
  addCcrMetric(ccrCoreExtractions, 1);

  const ui qSize = static_cast<ui>(Q.size());
  const ull universe = qSize == 64 ? ~0ULL
                                   : (qSize == 0 ? 0ULL
                                                 : (1ULL << qSize) - 1ULL);
  array<ull, 64> qAdj{};
  for (ui left = 0; left < qSize; ++left) {
    for (ui right = left + 1; right < qSize; ++right) {
      if (!adj(Q[left], Q[right]))
        continue;
      qAdj[left] |= 1ULL << right;
      qAdj[right] |= 1ULL << left;
    }
  }

  ull core = universe;
  ull residual = 0;
  ui coreCount = qSize;
  ui residualCount = 0;
  ull directedEdges = 0;
  for (ui local = 0; local < qSize; ++local)
    directedEdges += static_cast<ull>(
        __builtin_popcountll(qAdj[local] & core));

  while (coreCount > 1 &&
         directedEdges !=
             static_cast<ull>(coreCount) * (coreCount - 1)) {
    ui removed = numeric_limits<ui>::max();
    ui minimumDegree = numeric_limits<ui>::max();
    ull candidates = core;
    while (candidates != 0) {
      const ui local =
          static_cast<ui>(__builtin_ctzll(candidates));
      const ui degree = static_cast<ui>(
          __builtin_popcountll(qAdj[local] & core));
      if (degree < minimumDegree ||
          (degree == minimumDegree &&
           (removed == numeric_limits<ui>::max() ||
            Q[local] < Q[removed]))) {
        removed = local;
        minimumDegree = degree;
      }
      candidates &= candidates - 1;
    }
    directedEdges -= static_cast<ull>(minimumDegree) * 2;
    core &= ~(1ULL << removed);
    residual |= 1ULL << removed;
    --coreCount;
    ++residualCount;
  }
  addCcrMetric(ccrCoreVertices, coreCount);
  addCcrMetric(ccrResidualVertices, residualCount);

  const bool stopAfterOne = found != nullptr;
  auto recurse = [&](auto &&self, ull currentCore,
                     ull currentResidual, ull excluded,
                     ull selected, bool fromP) -> bool {
    if (stopAfterOne)
      addCcrMetric(ccrFindOneStates, 1);
    else
      addCcrMetric(ccrFullStates, 1);

    const size_t maximumSize =
        M.size() + static_cast<size_t>(__builtin_popcountll(
                       selected | currentCore | currentResidual));
    if (maximumSize < minCliqueSize)
      return false;

    // C vertices adjacent to all of P and P vertices adjacent to all of
    // C union P occur in every maximal extension. Move them directly to R.
    // Filtering X by the forced vertices preserves the usual BK invariant.
    ull forcedCore = 0;
    ull candidates = currentCore;
    while (candidates != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(candidates));
      const ull bit = 1ULL << local;
      if ((currentResidual & ~qAdj[local] & universe) == 0)
        forcedCore |= bit;
      candidates &= candidates - 1;
    }
    ull forcedResidual = 0;
    const ull active = currentCore | currentResidual;
    candidates = currentResidual;
    while (candidates != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(candidates));
      const ull bit = 1ULL << local;
      if (((active & ~qAdj[local] & universe) & ~bit) == 0)
        forcedResidual |= bit;
      candidates &= candidates - 1;
    }
    ull forced = forcedCore | forcedResidual;
    if (forced != 0) {
      currentCore &= ~forcedCore;
      currentResidual &= ~forcedResidual;
      selected |= forced;
      if (forcedResidual != 0)
        fromP = true;
      while (forced != 0) {
        const ui local = static_cast<ui>(__builtin_ctzll(forced));
        excluded &= qAdj[local];
        forced &= forced - 1;
      }
    }

    bool maximal = fromP &&
                   M.size() + static_cast<size_t>(__builtin_popcountll(
                                  selected | currentCore)) >=
                       minCliqueSize;
    auto localExtendsCore = [&](ui local) {
      return (currentCore & ~qAdj[local]) == 0;
    };
    if (maximal) {
      ull blockers = currentResidual | excluded;
      while (blockers != 0) {
        const ui local =
            static_cast<ui>(__builtin_ctzll(blockers));
        if (localExtendsCore(local)) {
          maximal = false;
          break;
        }
        blockers &= blockers - 1;
      }
    }
    if (maximal) {
      const ull represented = selected | currentCore;
      for (ui x : common) {
        bool extends = true;
        ull required = represented;
        while (required != 0) {
          const ui local =
              static_cast<ui>(__builtin_ctzll(required));
          if (!adj(x, Q[local])) {
            extends = false;
            break;
          }
          required &= required - 1;
        }
        if (extends) {
          maximal = false;
          break;
        }
      }
    }

    if (maximal) {
      ccrCliqueScratch.assign(M.begin(), M.end());
      ull extension = selected | currentCore;
      while (extension != 0) {
        const ui local =
            static_cast<ui>(__builtin_ctzll(extension));
        ccrCliqueScratch.push_back(Q[local]);
        extension &= extension - 1;
      }
      if (stopAfterOne) {
        *found = ccrCliqueScratch;
        sort(found->begin(), found->end());
        return true;
      }
      recordPureClique(ccrCliqueScratch);
    }
    if (currentResidual == 0)
      return false;

    // An excluded vertex with no P neighbor cannot survive any descendant,
    // because every descendant adds a vertex from P.
    ull excludedToCheck = excluded;
    while (excludedToCheck != 0) {
      const ui local =
          static_cast<ui>(__builtin_ctzll(excludedToCheck));
      const ull bit = 1ULL << local;
      if ((qAdj[local] & currentResidual) == 0)
        excluded &= ~bit;
      excludedToCheck &= excludedToCheck - 1;
    }

    // The sole P vertex is the only possible next non-core choice. Avoid
    // pivot selection and branch-root construction for this common terminal.
    if ((currentResidual & (currentResidual - 1)) == 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(currentResidual));
      const ull bit = 1ULL << local;
      return self(self, currentCore & qAdj[local], 0,
                  excluded & qAdj[local], selected | bit, true);
    }

    const ull pivotActive = currentCore | currentResidual;
    bool havePivot = false;
    ull pivotMask = 0;
    ui pivotLabel = 0;
    ui bestScore = 0;
    auto consider = [&](ull mask, ui label) {
      const ui score =
          static_cast<ui>(__builtin_popcountll(mask & pivotActive));
      if (!havePivot || score > bestScore ||
          (score == bestScore && label < pivotLabel)) {
        havePivot = true;
        pivotMask = mask;
        pivotLabel = label;
        bestScore = score;
      }
    };

    ull localPivots = pivotActive | excluded;
    while (localPivots != 0) {
      const ui local =
          static_cast<ui>(__builtin_ctzll(localPivots));
      consider(qAdj[local], Q[local]);
      localPivots &= localPivots - 1;
    }
    if (!havePivot)
      throw logic_error("small CCRMCE could not choose a pivot");

    ull roots = currentResidual & ~pivotMask & universe;
    ull coreNonNeighbors =
        currentCore & ~pivotMask & universe;
    while (coreNonNeighbors != 0) {
      const ui local =
          static_cast<ui>(__builtin_ctzll(coreNonNeighbors));
      if ((qAdj[local] & currentResidual & pivotMask) != 0)
        roots |= 1ULL << local;
      coreNonNeighbors &= coreNonNeighbors - 1;
    }

    while (roots != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(roots));
      const ull bit = 1ULL << local;
      const bool selectedFromP = (currentResidual & bit) != 0;
      if (self(self, currentCore & qAdj[local],
               currentResidual & qAdj[local],
               excluded & qAdj[local], selected | bit,
               selectedFromP) &&
          stopAfterOne)
        return true;

      if (selectedFromP)
        currentResidual &= ~bit;
      else
        currentCore &= ~bit;
      excluded |= bit;
      roots &= roots - 1;
    }
    return false;
  };

  return recurse(recurse, core, residual, 0, 0, true);
}

// Construct one reusable local induced graph for a formal branch (M,Q).
// Adjacency rows are dense only across Q; external X stays a compact list.
// This is the important implementation difference from the earlier
// correctness-first port: recursive CCRMCE states now use word operations
// instead of rebuilding sorted C/P/X vectors and issuing graph lookups.
void ReorderSib::prepareCcrBranch(const vector<ui> &Q,
                                  const vector<ui> &X,
                                  bool materializeExternalRows) {
  addCcrMetric(ccrCoreExtractions, 1);
  ccrQVertices.assign(Q.begin(), Q.end());
  ccrXVertices.assign(X.begin(), X.end());
  const size_t qSize = ccrQVertices.size();
  const size_t xSize = ccrXVertices.size();
  ccrWordCount = (qSize + 63) / 64;
  ccrExternalRowsMaterialized = materializeExternalRows;

#ifndef NDEBUG
  if (!is_sorted(Q.begin(), Q.end()) ||
      adjacent_find(Q.begin(), Q.end()) != Q.end() ||
      !is_sorted(X.begin(), X.end()) ||
      adjacent_find(X.begin(), X.end()) != X.end())
    throw logic_error("CCRMCE branch sets must be sorted and unique");
#endif

  if (++ccrIndexToken == 0) {
    fill(ccrIndexStamp.begin(), ccrIndexStamp.end(), 0);
    ccrIndexToken = 1;
  }
  for (size_t local = 0; local < qSize; ++local) {
    const ui vertex = ccrQVertices[local];
    ccrIndex[vertex] = static_cast<ui>(local);
    ccrIndexStamp[vertex] = ccrIndexToken;
  }

  if (ccrExternalRowsMaterialized) {
    for (size_t local = 0; local < xSize; ++local) {
      const ui vertex = ccrXVertices[local];
      ccrIndex[vertex] = static_cast<ui>(qSize + local);
      ccrIndexStamp[vertex] = ccrIndexToken;
    }
  }

  const size_t rowCount =
      qSize + (ccrExternalRowsMaterialized ? xSize : 0);
  if (ccrWordCount != 0 &&
      rowCount > numeric_limits<size_t>::max() / ccrWordCount)
    throw length_error("CCRMCE local adjacency dimensions overflow");
  ccrNeighborBits.assign(rowCount * ccrWordCount, 0);

  constexpr ull kAdjacencyLookupWeight = 8;
  ull scanWork = 0;
  for (ui vertex : ccrQVertices)
    scanWork += static_cast<ull>(adjacentVertices(vertex).size());
  ull pairWork = static_cast<ull>(qSize) *
                 static_cast<ull>(qSize > 0 ? qSize - 1 : 0) / 2;
  if (ccrExternalRowsMaterialized)
    pairWork += static_cast<ull>(qSize) * xSize;
  const ull globalLookupEquivalentScanWork =
      scanWork / kAdjacencyLookupWeight +
      static_cast<ull>(scanWork % kAdjacencyLookupWeight != 0);
  const bool scanRows =
      pairWork == 0 || globalLookupEquivalentScanWork <= pairWork;

  if (scanRows) {
    for (size_t local = 0; local < qSize; ++local) {
      ull *qRow = ccrNeighborBits.data() + local * ccrWordCount;
      for (ui neighbor : adjacentVertices(ccrQVertices[local])) {
        if (ccrIndexStamp[neighbor] == ccrIndexToken) {
          const size_t combinedLocal = ccrIndex[neighbor];
          if (combinedLocal < qSize) {
            qRow[combinedLocal >> 6] |=
                1ULL << (combinedLocal & 63);
            continue;
          }
          const size_t xLocal = combinedLocal - qSize;
          ull *xRow = ccrNeighborBits.data() +
                      (qSize + xLocal) * ccrWordCount;
          xRow[local >> 6] |= 1ULL << (local & 63);
        }
      }
    }
  } else {
    constexpr ull kPerRowAdjacencyLookupWeight = 6;
    ccrHeapPosition.assign(qSize, 0);
    for (size_t left = 0; left < qSize; ++left) {
      const AdjacencyRow row = adjacentVertices(ccrQVertices[left]);
      const ull lookupWork =
          static_cast<ull>(qSize - left - 1) +
          (ccrExternalRowsMaterialized ? static_cast<ull>(xSize) : 0ULL);
      const ull lookupEquivalentScanWork =
          static_cast<ull>(row.size()) / kPerRowAdjacencyLookupWeight +
          static_cast<ull>(row.size() % kPerRowAdjacencyLookupWeight != 0);
      const bool scanRow = lookupWork != 0 &&
                           lookupEquivalentScanWork <= lookupWork;
      ccrHeapPosition[left] = scanRow;

      if (scanRow) {
        for (ui neighbor : row) {
          if (ccrIndexStamp[neighbor] == ccrIndexToken) {
            const size_t combinedLocal = ccrIndex[neighbor];
            if (combinedLocal >= qSize) {
              const size_t xLocal = combinedLocal - qSize;
              ull *xRow = ccrNeighborBits.data() +
                          (qSize + xLocal) * ccrWordCount;
              xRow[left >> 6] |= 1ULL << (left & 63);
              continue;
            }
            const size_t right = combinedLocal;
            if (right <= left)
              continue;
            ull *leftRow =
                ccrNeighborBits.data() + left * ccrWordCount;
            ull *rightRow =
                ccrNeighborBits.data() + right * ccrWordCount;
            leftRow[right >> 6] |= 1ULL << (right & 63);
            rightRow[left >> 6] |= 1ULL << (left & 63);
          }
        }
        continue;
      }

      for (size_t right = left + 1; right < qSize; ++right) {
        if (!adj(ccrQVertices[left], ccrQVertices[right]))
          continue;
        ull *leftRow =
            ccrNeighborBits.data() + left * ccrWordCount;
        ull *rightRow =
            ccrNeighborBits.data() + right * ccrWordCount;
        leftRow[right >> 6] |= 1ULL << (right & 63);
        rightRow[left >> 6] |= 1ULL << (left & 63);
      }
    }

    if (ccrExternalRowsMaterialized) {
      for (size_t xLocal = 0; xLocal < xSize; ++xLocal) {
        ull *xRow = ccrNeighborBits.data() +
                    (qSize + xLocal) * ccrWordCount;
        for (size_t qLocal = 0; qLocal < qSize; ++qLocal) {
          if (ccrHeapPosition[qLocal] != 0)
            continue;
          if (adj(ccrXVertices[xLocal], ccrQVertices[qLocal]))
            xRow[qLocal >> 6] |= 1ULL << (qLocal & 63);
        }
      }
    }
  }

  if (ccrStates.size() < qSize + 1)
    ccrStates.resize(qSize + 1);
  CcrBitState &root = ccrStates[0];
  root.core.assign(ccrWordCount, ~0ULL);
  if (ccrWordCount != 0 && (qSize & 63) != 0)
    root.core.back() = (1ULL << (qSize & 63)) - 1;
  root.residual.assign(ccrWordCount, 0);
  root.excludedCandidates.assign(ccrWordCount, 0);
  root.excludedExternal.resize(xSize);
  iota(root.excludedExternal.begin(), root.excludedExternal.end(), 0);
  root.branchRoots.clear();
  root.coreCount = static_cast<ui>(qSize);
  root.residualCount = 0;

  if (qSize == 0) {
    addCcrMetric(ccrCoreVertices, 0);
    addCcrMetric(ccrResidualVertices, 0);
    return;
  }

  // A one-word induced graph is faster to peel by scanning its active bits
  // than by maintaining three heap vectors. This chooses the same
  // minimum-degree vertex (including the original vertex-ID tie break), so
  // only branch-preparation cost changes; the recursive CCRMCE state and
  // external-X handling remain identical.
  if (ccrWordCount == 1) {
    ull core = root.core[0];
    ull directedEdges = 0;
    for (ui local = 0; local < qSize; ++local)
      directedEdges += static_cast<ull>(__builtin_popcountll(
          ccrNeighborBits[local] & core));

    while (root.coreCount > 1 &&
           directedEdges !=
               static_cast<ull>(root.coreCount) * (root.coreCount - 1)) {
      ui removed = numeric_limits<ui>::max();
      ui minimumDegree = numeric_limits<ui>::max();
      ull candidates = core;
      while (candidates != 0) {
        const ui local = static_cast<ui>(__builtin_ctzll(candidates));
        const ui degree = static_cast<ui>(__builtin_popcountll(
            ccrNeighborBits[local] & core));
        if (degree < minimumDegree ||
            (degree == minimumDegree &&
             (removed == numeric_limits<ui>::max() ||
              ccrQVertices[local] < ccrQVertices[removed]))) {
          removed = local;
          minimumDegree = degree;
        }
        candidates &= candidates - 1;
      }
      directedEdges -= static_cast<ull>(minimumDegree) * 2;
      const ull bit = 1ULL << removed;
      core &= ~bit;
      root.core[0] = core;
      root.residual[0] |= bit;
      --root.coreCount;
      ++root.residualCount;
    }

    addCcrMetric(ccrCoreVertices, root.coreCount);
    addCcrMetric(ccrResidualVertices, root.residualCount);
    return;
  }

  ccrDegree.resize(qSize);
  ccrHeap.resize(qSize);
  ccrHeapPosition.resize(qSize);
  ull directedEdges = 0;
  for (ui local = 0; local < qSize; ++local) {
    const ull *row =
        ccrNeighborBits.data() + static_cast<size_t>(local) * ccrWordCount;
    ui degree = 0;
    for (size_t word = 0; word < ccrWordCount; ++word)
      degree += static_cast<ui>(__builtin_popcountll(row[word]));
    ccrDegree[local] = degree;
    ccrHeap[local] = local;
    ccrHeapPosition[local] = local;
    directedEdges += degree;
  }

  auto heapLess = [&](ui lhs, ui rhs) {
    if (ccrDegree[lhs] != ccrDegree[rhs])
      return ccrDegree[lhs] < ccrDegree[rhs];
    return ccrQVertices[lhs] < ccrQVertices[rhs];
  };
  auto heapSwap = [&](size_t lhs, size_t rhs) {
    swap(ccrHeap[lhs], ccrHeap[rhs]);
    ccrHeapPosition[ccrHeap[lhs]] = static_cast<ui>(lhs);
    ccrHeapPosition[ccrHeap[rhs]] = static_cast<ui>(rhs);
  };
  auto siftDown = [&](size_t position) {
    while (true) {
      const size_t left = position * 2 + 1;
      if (left >= ccrHeap.size())
        break;
      size_t best = left;
      const size_t right = left + 1;
      if (right < ccrHeap.size() &&
          heapLess(ccrHeap[right], ccrHeap[left]))
        best = right;
      if (!heapLess(ccrHeap[best], ccrHeap[position]))
        break;
      heapSwap(position, best);
      position = best;
    }
  };
  auto siftUp = [&](size_t position) {
    while (position != 0) {
      const size_t parent = (position - 1) / 2;
      if (!heapLess(ccrHeap[position], ccrHeap[parent]))
        break;
      heapSwap(position, parent);
      position = parent;
    }
  };
  for (size_t position = ccrHeap.size() / 2; position-- > 0;)
    siftDown(position);

  size_t activeCount = qSize;
  while (activeCount > 1) {
    ull completeDirectedEdges = 0;
    if (!tryMultiplyUll(static_cast<ull>(activeCount),
                        static_cast<ull>(activeCount - 1),
                        completeDirectedEdges))
      throw overflow_error("CCRMCE core completeness test exceeds uint64_t");
    if (directedEdges == completeDirectedEdges)
      break;

    const ui removed = ccrHeap.front();
    heapSwap(0, ccrHeap.size() - 1);
    ccrHeap.pop_back();
    ccrHeapPosition[removed] = numeric_limits<ui>::max();
    if (!ccrHeap.empty())
      siftDown(0);

    const ull removedDirectedEdges =
        static_cast<ull>(ccrDegree[removed]) * 2;
    if (removedDirectedEdges > directedEdges)
      throw logic_error("CCRMCE induced edge count became inconsistent");
    directedEdges -= removedDirectedEdges;
    ccrClearBit(root.core, removed);
    ccrSetBit(root.residual, removed);
    --root.coreCount;
    ++root.residualCount;
    --activeCount;

    const ull *row = ccrNeighborBits.data() +
                     static_cast<size_t>(removed) * ccrWordCount;
    for (size_t wordIndex = 0; wordIndex < ccrWordCount; ++wordIndex) {
      ull neighbors = row[wordIndex] & root.core[wordIndex];
      while (neighbors != 0) {
        const unsigned bit =
            static_cast<unsigned>(__builtin_ctzll(neighbors));
        const ui neighbor =
            static_cast<ui>(wordIndex * 64 + bit);
        if (ccrDegree[neighbor] == 0)
          throw logic_error("CCRMCE active-neighbor degree underflow");
        --ccrDegree[neighbor];
        siftUp(ccrHeapPosition[neighbor]);
        neighbors &= neighbors - 1;
      }
    }
  }

  addCcrMetric(ccrCoreVertices, root.coreCount);
  addCcrMetric(ccrResidualVertices, root.residualCount);
}

void ReorderSib::normalizeCcrState(vector<ui> &R, CcrBitState &state,
                                   bool &fromP) {
  state.branchRoots.clear();
  state.branchRoots.reserve(
      static_cast<size_t>(state.coreCount) + state.residualCount);

  // A core vertex adjacent to every residual candidate belongs to every
  // maximal extension of this state: C is already a clique, so omitting the
  // vertex would leave the extension non-maximal.
  ccrForEachBit(state.core, [&](ui local) {
    const ull *row = ccrNeighborBits.data() +
                     static_cast<size_t>(local) * ccrWordCount;
    bool universal = true;
    for (size_t word = 0; word < ccrWordCount; ++word) {
      if ((state.residual[word] & ~row[word]) != 0) {
        universal = false;
        break;
      }
    }
    if (universal)
      state.branchRoots.push_back(local);
  });
  const size_t forcedCoreCount = state.branchRoots.size();

  // A residual vertex adjacent to all other active vertices is likewise
  // forced. The vertex's own bit is ignored because adjacency has no loops.
  ccrForEachBit(state.residual, [&](ui local) {
    const ull *row = ccrNeighborBits.data() +
                     static_cast<size_t>(local) * ccrWordCount;
    bool universal = true;
    for (size_t word = 0; word < ccrWordCount; ++word) {
      ull nonNeighbors =
          (state.core[word] | state.residual[word]) & ~row[word];
      if (word == (local >> 6))
        nonNeighbors &= ~(1ULL << (local & 63));
      if (nonNeighbors != 0) {
        universal = false;
        break;
      }
    }
    if (universal)
      state.branchRoots.push_back(local);
  });

  if (state.branchRoots.empty())
    return;

  for (size_t index = 0; index < state.branchRoots.size(); ++index) {
    const ui local = state.branchRoots[index];
    if (index < forcedCoreCount) {
      ccrClearBit(state.core, local);
      --state.coreCount;
    } else {
      ccrClearBit(state.residual, local);
      --state.residualCount;
      fromP = true;
    }
    R.push_back(ccrQVertices[local]);

    const ull *row = ccrNeighborBits.data() +
                     static_cast<size_t>(local) * ccrWordCount;
    for (size_t word = 0; word < ccrWordCount; ++word)
      state.excludedCandidates[word] &= row[word];
  }

  size_t kept = 0;
  const size_t qSize = ccrQVertices.size();
  for (ui xLocal : state.excludedExternal) {
    bool adjacentToAll = true;
    if (ccrExternalRowsMaterialized) {
      const ull *xRow = ccrNeighborBits.data() +
                        (qSize + xLocal) * ccrWordCount;
      for (ui local : state.branchRoots) {
        if ((xRow[local >> 6] & (1ULL << (local & 63))) == 0) {
          adjacentToAll = false;
          break;
        }
      }
    } else {
      for (ui local : state.branchRoots) {
        if (!adj(ccrXVertices[xLocal], ccrQVertices[local])) {
          adjacentToAll = false;
          break;
        }
      }
    }
    if (adjacentToAll)
      state.excludedExternal[kept++] = xLocal;
  }
  state.excludedExternal.resize(kept);
  state.branchRoots.clear();
}

void ReorderSib::pruneCcrExcludedWithoutResidualNeighbors(
    CcrBitState &state) {
  ccrForEachBit(state.excludedCandidates, [&](ui local) {
    const ull *row = ccrNeighborBits.data() +
                     static_cast<size_t>(local) * ccrWordCount;
    bool hasResidualNeighbor = false;
    for (size_t word = 0; word < ccrWordCount; ++word) {
      if ((row[word] & state.residual[word]) != 0) {
        hasResidualNeighbor = true;
        break;
      }
    }
    if (!hasResidualNeighbor)
      ccrClearBit(state.excludedCandidates, local);
  });

  size_t kept = 0;
  const size_t qSize = ccrQVertices.size();
  for (ui xLocal : state.excludedExternal) {
    bool hasResidualNeighbor = false;
    if (ccrExternalRowsMaterialized) {
      const ull *xRow = ccrNeighborBits.data() +
                        (qSize + xLocal) * ccrWordCount;
      for (size_t word = 0; word < ccrWordCount; ++word) {
        if ((xRow[word] & state.residual[word]) != 0) {
          hasResidualNeighbor = true;
          break;
        }
      }
    } else {
      for (size_t wordIndex = 0;
           wordIndex < ccrWordCount && !hasResidualNeighbor; ++wordIndex) {
        ull residual = state.residual[wordIndex];
        while (residual != 0) {
          const ui local = static_cast<ui>(
              wordIndex * 64 + __builtin_ctzll(residual));
          if (adj(ccrXVertices[xLocal], ccrQVertices[local])) {
            hasResidualNeighbor = true;
            break;
          }
          residual &= residual - 1;
        }
      }
    }
    if (hasResidualNeighbor)
      state.excludedExternal[kept++] = xLocal;
  }
  state.excludedExternal.resize(kept);
}

void ReorderSib::buildCcrChildState(CcrBitState &child,
                                    const CcrBitState &parent,
                                    ui local) const {
  child.core.resize(ccrWordCount);
  child.residual.resize(ccrWordCount);
  child.excludedCandidates.resize(ccrWordCount);
  const ull *row = ccrNeighborBits.data() +
                   static_cast<size_t>(local) * ccrWordCount;
  for (size_t word = 0; word < ccrWordCount; ++word) {
    child.core[word] = parent.core[word] & row[word];
    child.residual[word] = parent.residual[word] & row[word];
    child.excludedCandidates[word] =
        parent.excludedCandidates[word] & row[word];
  }
  child.coreCount = ccrBitCount(child.core);
  child.residualCount = ccrBitCount(child.residual);
  child.excludedExternal.clear();
  if (child.excludedExternal.capacity() <
      parent.excludedExternal.size())
    child.excludedExternal.reserve(parent.excludedExternal.size());
  const size_t qSize = ccrQVertices.size();
  for (ui xLocal : parent.excludedExternal) {
    const bool adjacent =
        ccrExternalRowsMaterialized
            ? (ccrNeighborBits[(qSize + xLocal) * ccrWordCount +
                               (local >> 6)] &
               (1ULL << (local & 63))) != 0
            : adj(ccrXVertices[xLocal], ccrQVertices[local]);
    if (adjacent)
      child.excludedExternal.push_back(xLocal);
  }
  child.branchRoots.clear();
}

bool ReorderSib::ccrCoreUnionIsMaximal(
    const CcrBitState &state) const {
  auto rowContainsCore = [&](size_t rowIndex) {
    const ull *row =
        ccrNeighborBits.data() + rowIndex * ccrWordCount;
    for (size_t word = 0; word < ccrWordCount; ++word)
      if ((state.core[word] & ~row[word]) != 0)
        return false;
    return true;
  };

  bool blocked = false;
  ccrForEachBit(state.residual, [&](ui local) {
    blocked = blocked || rowContainsCore(local);
  });
  if (blocked)
    return false;
  ccrForEachBit(state.excludedCandidates, [&](ui local) {
    blocked = blocked || rowContainsCore(local);
  });
  if (blocked)
    return false;
  const size_t qSize = ccrQVertices.size();
  for (ui xLocal : state.excludedExternal) {
    if (ccrExternalRowsMaterialized) {
      if (rowContainsCore(qSize + xLocal))
        return false;
      continue;
    }
    bool extendsCore = true;
    ccrForEachBit(state.core, [&](ui local) {
      extendsCore =
          extendsCore &&
          adj(ccrXVertices[xLocal], ccrQVertices[local]);
    });
    if (extendsCore)
      return false;
  }
  return true;
}

size_t ReorderSib::selectCcrPivot(const CcrBitState &state) const {
  const size_t qSize = ccrQVertices.size();
  bool havePivot = false;
  size_t pivot = 0;
  ui pivotLabel = 0;
  ui bestScore = 0;

  auto consider = [&](size_t rowIndex, ui label) {
    const ull *row =
        ccrNeighborBits.data() + rowIndex * ccrWordCount;
    ui score = 0;
    for (size_t word = 0; word < ccrWordCount; ++word) {
      score += static_cast<ui>(__builtin_popcountll(
          row[word] & (state.core[word] | state.residual[word])));
    }
    if (!havePivot || score > bestScore ||
        (score == bestScore && label < pivotLabel)) {
      havePivot = true;
      pivot = rowIndex;
      pivotLabel = label;
      bestScore = score;
    }
  };

  for (size_t wordIndex = 0; wordIndex < ccrWordCount; ++wordIndex) {
    ull localCandidates = state.core[wordIndex] |
                          state.residual[wordIndex] |
                          state.excludedCandidates[wordIndex];
    while (localCandidates != 0) {
      const unsigned bit =
          static_cast<unsigned>(__builtin_ctzll(localCandidates));
      const ui local = static_cast<ui>(wordIndex * 64 + bit);
      consider(local, ccrQVertices[local]);
      localCandidates &= localCandidates - 1;
    }
  }
  if (ccrExternalRowsMaterialized)
    for (ui xLocal : state.excludedExternal)
      consider(qSize + xLocal, ccrXVertices[xLocal]);

  if (!havePivot)
    throw logic_error("CCRMCE pivot requested for an empty state");
  return pivot;
}

void ReorderSib::buildCcrBranchRoots(CcrBitState &state,
                                     size_t pivot) const {
  state.branchRoots.clear();
  state.branchRoots.reserve(
      static_cast<size_t>(state.coreCount) + state.residualCount);
  const ull *pivotRow =
      ccrNeighborBits.data() + pivot * ccrWordCount;

  for (size_t wordIndex = 0; wordIndex < ccrWordCount; ++wordIndex) {
    ull roots = state.residual[wordIndex] & ~pivotRow[wordIndex];
    while (roots != 0) {
      const unsigned bit = static_cast<unsigned>(__builtin_ctzll(roots));
      state.branchRoots.push_back(
          static_cast<ui>(wordIndex * 64 + bit));
      roots &= roots - 1;
    }
  }

  for (size_t wordIndex = 0; wordIndex < ccrWordCount; ++wordIndex) {
    ull roots = state.core[wordIndex] & ~pivotRow[wordIndex];
    while (roots != 0) {
      const unsigned bit = static_cast<unsigned>(__builtin_ctzll(roots));
      const ui local = static_cast<ui>(wordIndex * 64 + bit);
      const ull *row = ccrNeighborBits.data() +
                       static_cast<size_t>(local) * ccrWordCount;
      bool sharesPivotResidualNeighbor = false;
      for (size_t word = 0; word < ccrWordCount; ++word) {
        if ((row[word] & state.residual[word] & pivotRow[word]) != 0) {
          sharesPivotResidualNeighbor = true;
          break;
        }
      }
      if (sharesPivotResidualNeighbor)
        state.branchRoots.push_back(local);
      roots &= roots - 1;
    }
  }
}

// Full CCRMCE recursion for a materialized one-word branch. prepareCcrBranch
// has already built every Q and external-X row, so keeping C/P/X-within-Q in
// scalar masks removes per-state dynamic bit-vector work without returning to
// the slower lazy-X scalar path.
void ReorderSib::enumeratePreparedOneWordBranch(const vector<ui> &M) {
  if (ccrWordCount != 1 || !ccrExternalRowsMaterialized)
    throw logic_error("prepared one-word CCRMCE requires materialized rows");

  const size_t qSize = ccrQVertices.size();
  const ull universe =
      qSize == 64 ? ~0ULL : ((1ULL << static_cast<unsigned>(qSize)) - 1ULL);
  auto qRow = [&](ui local) {
    return ccrNeighborBits[static_cast<size_t>(local)];
  };
  auto xRow = [&](ui xLocal) {
    return ccrNeighborBits[qSize + static_cast<size_t>(xLocal)];
  };

  auto recurse = [&](auto &&self, ull core, ull residual,
                     ull excluded, ull selected, bool fromP,
                     vector<ui> &external, size_t depth) -> void {
    addCcrMetric(ccrFullStates, 1);
    if (M.size() + static_cast<size_t>(__builtin_popcountll(
                       selected | core | residual)) <
        minCliqueSize)
      return;

    ull forcedCore = 0;
    ull candidates = core;
    while (candidates != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(candidates));
      const ull bit = 1ULL << local;
      if ((residual & ~qRow(local) & universe) == 0)
        forcedCore |= bit;
      candidates &= candidates - 1;
    }

    ull forcedResidual = 0;
    const ull active = core | residual;
    candidates = residual;
    while (candidates != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(candidates));
      const ull bit = 1ULL << local;
      if (((active & ~qRow(local) & universe) & ~bit) == 0)
        forcedResidual |= bit;
      candidates &= candidates - 1;
    }

    const ull forced = forcedCore | forcedResidual;
    if (forced != 0) {
      core &= ~forcedCore;
      residual &= ~forcedResidual;
      selected |= forced;
      if (forcedResidual != 0)
        fromP = true;
      ull forcedBits = forced;
      while (forcedBits != 0) {
        const ui local = static_cast<ui>(__builtin_ctzll(forcedBits));
        excluded &= qRow(local);
        forcedBits &= forcedBits - 1;
      }

      size_t kept = 0;
      for (ui xLocal : external)
        if ((xRow(xLocal) & forced) == forced)
          external[kept++] = xLocal;
      external.resize(kept);
    }

    bool maximal =
        fromP &&
        M.size() + static_cast<size_t>(
                       __builtin_popcountll(selected | core)) >=
            minCliqueSize;
    if (maximal) {
      ull blockers = residual | excluded;
      while (blockers != 0) {
        const ui local = static_cast<ui>(__builtin_ctzll(blockers));
        if ((core & ~qRow(local)) == 0) {
          maximal = false;
          break;
        }
        blockers &= blockers - 1;
      }
    }
    if (maximal)
      for (ui xLocal : external)
        if ((core & ~xRow(xLocal)) == 0) {
          maximal = false;
          break;
        }

    if (maximal) {
      ccrCliqueScratch.assign(M.begin(), M.end());
      ull extension = selected | core;
      while (extension != 0) {
        const ui local = static_cast<ui>(__builtin_ctzll(extension));
        ccrCliqueScratch.push_back(ccrQVertices[local]);
        extension &= extension - 1;
      }
      recordPureClique(ccrCliqueScratch);
    }
    if (residual == 0)
      return;

    ull excludedToCheck = excluded;
    while (excludedToCheck != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(excludedToCheck));
      const ull bit = 1ULL << local;
      if ((qRow(local) & residual) == 0)
        excluded &= ~bit;
      excludedToCheck &= excludedToCheck - 1;
    }
    size_t kept = 0;
    for (ui xLocal : external)
      if ((xRow(xLocal) & residual) != 0)
        external[kept++] = xLocal;
    external.resize(kept);

    auto descend = [&](ui local, bool selectedFromP) {
      vector<ui> &childExternal = ccrStates[depth + 1].excludedExternal;
      childExternal.clear();
      if (childExternal.capacity() < external.size())
        childExternal.reserve(external.size());
      const ull bit = 1ULL << local;
      for (ui xLocal : external)
        if ((xRow(xLocal) & bit) != 0)
          childExternal.push_back(xLocal);
      const ull row = qRow(local);
      self(self, core & row, residual & row, excluded & row,
           selected | bit, selectedFromP, childExternal, depth + 1);
    };

    if ((residual & (residual - 1)) == 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(residual));
      descend(local, true);
      return;
    }

    const ull pivotActive = core | residual;
    bool havePivot = false;
    ull pivotMask = 0;
    ui pivotLabel = 0;
    ui bestScore = 0;
    auto consider = [&](ull row, ui label) {
      const ui score =
          static_cast<ui>(__builtin_popcountll(row & pivotActive));
      if (!havePivot || score > bestScore ||
          (score == bestScore && label < pivotLabel)) {
        havePivot = true;
        pivotMask = row;
        pivotLabel = label;
        bestScore = score;
      }
    };

    ull localPivots = pivotActive | excluded;
    while (localPivots != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(localPivots));
      consider(qRow(local), ccrQVertices[local]);
      localPivots &= localPivots - 1;
    }
    for (ui xLocal : external)
      consider(xRow(xLocal), ccrXVertices[xLocal]);
    if (!havePivot)
      throw logic_error("prepared one-word CCRMCE could not choose a pivot");

    ull roots = residual & ~pivotMask & universe;
    ull coreNonNeighbors = core & ~pivotMask & universe;
    while (coreNonNeighbors != 0) {
      const ui local =
          static_cast<ui>(__builtin_ctzll(coreNonNeighbors));
      if ((qRow(local) & residual & pivotMask) != 0)
        roots |= 1ULL << local;
      coreNonNeighbors &= coreNonNeighbors - 1;
    }

    while (roots != 0) {
      const ui local = static_cast<ui>(__builtin_ctzll(roots));
      const ull bit = 1ULL << local;
      const bool selectedFromP = (residual & bit) != 0;
      descend(local, selectedFromP);
      if (selectedFromP)
        residual &= ~bit;
      else
        core &= ~bit;
      excluded |= bit;
      roots &= roots - 1;
    }
  };

  CcrBitState &root = ccrStates[0];
  recurse(recurse, root.core[0], root.residual[0],
          root.excludedCandidates[0], 0, true,
          root.excludedExternal, 0);
}

// The same materialized scalar-state specialization for 65--128 candidates.
// Two fixed words avoid allocating and copying three vector objects at each
// recursive state while retaining the exact CCRMCE reductions and T1/T2 set.
void ReorderSib::enumeratePreparedTwoWordBranch(const vector<ui> &M) {
  if (ccrWordCount != 2 || !ccrExternalRowsMaterialized)
    throw logic_error("prepared two-word CCRMCE requires materialized rows");

  using Mask = array<ull, 2>;
  const size_t qSize = ccrQVertices.size();
  const Mask universe{
      ~0ULL,
      (qSize & 63) == 0
          ? ~0ULL
          : (1ULL << static_cast<unsigned>(qSize & 63)) - 1ULL};
  auto row = [&](size_t rowIndex) {
    return ccrNeighborBits.data() + rowIndex * 2;
  };
  auto count = [](const Mask &mask) {
    return static_cast<ui>(__builtin_popcountll(mask[0]) +
                           __builtin_popcountll(mask[1]));
  };
  auto forEach = [](const Mask &mask, auto &&visit) {
    for (size_t word = 0; word < 2; ++word) {
      ull bits = mask[word];
      while (bits != 0) {
        const ui local = static_cast<ui>(
            word * 64 + static_cast<size_t>(__builtin_ctzll(bits)));
        visit(local);
        bits &= bits - 1;
      }
    }
  };

  auto recurse = [&](auto &&self, Mask core, Mask residual,
                     Mask excluded, Mask selected, bool fromP,
                     vector<ui> &external, size_t depth) -> void {
    addCcrMetric(ccrFullStates, 1);
    const Mask possible{selected[0] | core[0] | residual[0],
                        selected[1] | core[1] | residual[1]};
    if (M.size() + count(possible) < minCliqueSize)
      return;

    Mask forcedCore{};
    forEach(core, [&](ui local) {
      const ull *localRow = row(local);
      if ((residual[0] & ~localRow[0]) == 0 &&
          (residual[1] & ~localRow[1]) == 0)
        forcedCore[local >> 6] |= 1ULL << (local & 63);
    });

    const Mask active{core[0] | residual[0], core[1] | residual[1]};
    Mask forcedResidual{};
    forEach(residual, [&](ui local) {
      const ull *localRow = row(local);
      ull nonNeighbors0 = active[0] & ~localRow[0];
      ull nonNeighbors1 = active[1] & ~localRow[1];
      if ((local >> 6) == 0)
        nonNeighbors0 &= ~(1ULL << (local & 63));
      else
        nonNeighbors1 &= ~(1ULL << (local & 63));
      if (nonNeighbors0 == 0 && nonNeighbors1 == 0)
        forcedResidual[local >> 6] |= 1ULL << (local & 63);
    });

    const Mask forced{forcedCore[0] | forcedResidual[0],
                      forcedCore[1] | forcedResidual[1]};
    if ((forced[0] | forced[1]) != 0) {
      for (size_t word = 0; word < 2; ++word) {
        core[word] &= ~forcedCore[word];
        residual[word] &= ~forcedResidual[word];
        selected[word] |= forced[word];
      }
      if ((forcedResidual[0] | forcedResidual[1]) != 0)
        fromP = true;
      forEach(forced, [&](ui local) {
        const ull *localRow = row(local);
        excluded[0] &= localRow[0];
        excluded[1] &= localRow[1];
      });

      size_t kept = 0;
      for (ui xLocal : external) {
        const ull *externalRow = row(qSize + xLocal);
        if ((forced[0] & ~externalRow[0]) == 0 &&
            (forced[1] & ~externalRow[1]) == 0)
          external[kept++] = xLocal;
      }
      external.resize(kept);
    }

    const Mask represented{selected[0] | core[0],
                           selected[1] | core[1]};
    bool maximal =
        fromP && M.size() + count(represented) >= minCliqueSize;
    if (maximal) {
      const Mask blockers{residual[0] | excluded[0],
                          residual[1] | excluded[1]};
      forEach(blockers, [&](ui local) {
        if (!maximal)
          return;
        const ull *localRow = row(local);
        if ((core[0] & ~localRow[0]) == 0 &&
            (core[1] & ~localRow[1]) == 0)
          maximal = false;
      });
    }
    if (maximal)
      for (ui xLocal : external) {
        const ull *externalRow = row(qSize + xLocal);
        if ((core[0] & ~externalRow[0]) == 0 &&
            (core[1] & ~externalRow[1]) == 0) {
          maximal = false;
          break;
        }
      }

    if (maximal) {
      ccrCliqueScratch.assign(M.begin(), M.end());
      forEach(represented, [&](ui local) {
        ccrCliqueScratch.push_back(ccrQVertices[local]);
      });
      recordPureClique(ccrCliqueScratch);
    }
    if ((residual[0] | residual[1]) == 0)
      return;

    const Mask excludedToCheck = excluded;
    forEach(excludedToCheck, [&](ui local) {
      const ull *localRow = row(local);
      if (((localRow[0] & residual[0]) |
           (localRow[1] & residual[1])) == 0)
        excluded[local >> 6] &= ~(1ULL << (local & 63));
    });
    size_t kept = 0;
    for (ui xLocal : external) {
      const ull *externalRow = row(qSize + xLocal);
      if (((externalRow[0] & residual[0]) |
           (externalRow[1] & residual[1])) != 0)
        external[kept++] = xLocal;
    }
    external.resize(kept);

    auto descend = [&](ui local, bool selectedFromP) {
      vector<ui> &childExternal = ccrStates[depth + 1].excludedExternal;
      childExternal.clear();
      if (childExternal.capacity() < external.size())
        childExternal.reserve(external.size());
      const size_t localWord = local >> 6;
      const ull bit = 1ULL << (local & 63);
      for (ui xLocal : external)
        if ((row(qSize + xLocal)[localWord] & bit) != 0)
          childExternal.push_back(xLocal);
      const ull *localRow = row(local);
      const Mask childCore{core[0] & localRow[0],
                           core[1] & localRow[1]};
      const Mask childResidual{residual[0] & localRow[0],
                               residual[1] & localRow[1]};
      const Mask childExcluded{excluded[0] & localRow[0],
                               excluded[1] & localRow[1]};
      Mask childSelected = selected;
      childSelected[localWord] |= bit;
      self(self, childCore, childResidual, childExcluded,
           childSelected, selectedFromP, childExternal, depth + 1);
    };

    if (count(residual) == 1) {
      const ui local = residual[0] != 0
                           ? static_cast<ui>(__builtin_ctzll(residual[0]))
                           : static_cast<ui>(64 +
                                             __builtin_ctzll(residual[1]));
      descend(local, true);
      return;
    }

    const Mask pivotActive{core[0] | residual[0],
                           core[1] | residual[1]};
    bool havePivot = false;
    const ull *pivotRow = nullptr;
    ui pivotLabel = 0;
    ui bestScore = 0;
    auto consider = [&](const ull *candidateRow, ui label) {
      const ui score = static_cast<ui>(
          __builtin_popcountll(candidateRow[0] & pivotActive[0]) +
          __builtin_popcountll(candidateRow[1] & pivotActive[1]));
      if (!havePivot || score > bestScore ||
          (score == bestScore && label < pivotLabel)) {
        havePivot = true;
        pivotRow = candidateRow;
        pivotLabel = label;
        bestScore = score;
      }
    };

    const Mask localPivots{pivotActive[0] | excluded[0],
                           pivotActive[1] | excluded[1]};
    forEach(localPivots, [&](ui local) {
      consider(row(local), ccrQVertices[local]);
    });
    for (ui xLocal : external)
      consider(row(qSize + xLocal), ccrXVertices[xLocal]);
    if (!havePivot)
      throw logic_error("prepared two-word CCRMCE could not choose a pivot");

    Mask roots{residual[0] & ~pivotRow[0] & universe[0],
               residual[1] & ~pivotRow[1] & universe[1]};
    const Mask coreNonNeighbors{core[0] & ~pivotRow[0] & universe[0],
                                core[1] & ~pivotRow[1] & universe[1]};
    forEach(coreNonNeighbors, [&](ui local) {
      const ull *localRow = row(local);
      if (((localRow[0] & residual[0] & pivotRow[0]) |
           (localRow[1] & residual[1] & pivotRow[1])) != 0)
        roots[local >> 6] |= 1ULL << (local & 63);
    });

    for (size_t word = 0; word < 2; ++word) {
      while (roots[word] != 0) {
        const ui local = static_cast<ui>(
            word * 64 + static_cast<size_t>(__builtin_ctzll(roots[word])));
        const ull bit = 1ULL << (local & 63);
        const bool selectedFromP = (residual[word] & bit) != 0;
        descend(local, selectedFromP);
        if (selectedFromP)
          residual[word] &= ~bit;
        else
          core[word] &= ~bit;
        excluded[word] |= bit;
        roots[word] &= roots[word] - 1;
      }
    }
  };

  CcrBitState &root = ccrStates[0];
  const Mask rootCore{root.core[0], root.core[1]};
  const Mask rootResidual{root.residual[0], root.residual[1]};
  const Mask rootExcluded{root.excludedCandidates[0],
                          root.excludedCandidates[1]};
  const Mask rootSelected{};
  recurse(recurse, rootCore, rootResidual, rootExcluded,
          rootSelected, true, root.excludedExternal, 0);
}

bool ReorderSib::findOnePure(const vector<ui> &M, const vector<ui> &Q,
                             vector<ui> &found) {
  found.clear();
  if (M.empty())
    return false;

  const AdjacencyRow firstRow = adjacentVertices(M[0]);
  const size_t forwardOffset = firstForwardNeighbor[M[0]];
  const bool wholeRootBranch =
      M.size() == 1 && Q.size() == firstRow.size() - forwardOffset;
  const vector<ui> *external = nullptr;
  if (wholeRootBranch) {
    ccrBranchX.assign(firstRow.begin(), firstRow.begin() + forwardOffset);
    external = &ccrBranchX;
  } else {
    ccrCommon.assign(firstRow.begin(), firstRow.end());
    for (size_t i = 1; i < M.size() && !ccrCommon.empty(); ++i) {
      intersectInto(ccrCommonScratch, ccrCommon, adjacentVertices(M[i]));
      ccrCommon.swap(ccrCommonScratch);
    }
  }
  addCcrMetric(ccrFindOneCalls, 1);
  if (Q.size() <= 4)
    return runSmallCcrBranch(M, Q,
                             wholeRootBranch ? *external : ccrCommon,
                             &found);
  if (!wholeRootBranch)
    setDiffInto(ccrBranchX, ccrCommon, Q);
  prepareCcrBranch(Q, ccrBranchX, false);
  ccrBranchR.assign(M.begin(), M.end());
  return findOnePureRecursive(ccrBranchR, ccrStates[0], true, found, 0);
}

bool ReorderSib::findOnePureRecursive(vector<ui> &R,
                                      CcrBitState &state, bool fromP,
                                      vector<ui> &found, size_t depth) {
  addCcrMetric(ccrFindOneStates, 1);
  if (R.size() + state.coreCount + state.residualCount < minCliqueSize)
    return false;

  const size_t entryRSize = R.size();
  normalizeCcrState(R, state, fromP);

  if (fromP && R.size() + state.coreCount >= minCliqueSize &&
      ccrCoreUnionIsMaximal(state)) {
    found = R;
    ccrForEachBit(state.core, [&](ui local) {
      found.push_back(ccrQVertices[local]);
    });
    sort(found.begin(), found.end());
    R.resize(entryRSize);
    return true;
  }
  if (state.residualCount == 0) {
    R.resize(entryRSize);
    return false;
  }

  pruneCcrExcludedWithoutResidualNeighbors(state);
  if (state.residualCount == 1) {
    const ui local = ccrFirstBit(state.residual);
    CcrBitState &child = ccrStates[depth + 1];
    buildCcrChildState(child, state, local);
    R.push_back(ccrQVertices[local]);
    const bool result =
        findOnePureRecursive(R, child, true, found, depth + 1);
    R.resize(entryRSize);
    return result;
  }

  const size_t pivot = selectCcrPivot(state);
  buildCcrBranchRoots(state, pivot);
  for (ui local : state.branchRoots) {
    const bool selectedFromP = ccrBitIsSet(state.residual, local);
    if (!selectedFromP && !ccrBitIsSet(state.core, local))
      throw logic_error("CCRMCE find-one branch root was already removed");

    CcrBitState &child = ccrStates[depth + 1];
    buildCcrChildState(child, state, local);

    R.push_back(ccrQVertices[local]);
    if (findOnePureRecursive(R, child, selectedFromP, found,
                             depth + 1)) {
      R.resize(entryRSize);
      return true;
    }
    R.pop_back();

    if (selectedFromP) {
      ccrClearBit(state.residual, local);
      --state.residualCount;
    } else {
      ccrClearBit(state.core, local);
      --state.coreCount;
    }
    ccrSetBit(state.excludedCandidates, local);
  }
  R.resize(entryRSize);
  return false;
}

void ReorderSib::enumerateAllPureBranch(const vector<ui> &M,
                                        const vector<ui> &Q) {
  if (M.empty())
    return;

  const AdjacencyRow firstRow = adjacentVertices(M[0]);
  const size_t forwardOffset = firstForwardNeighbor[M[0]];
  const bool wholeRootBranch =
      M.size() == 1 && Q.size() == firstRow.size() - forwardOffset;
  const vector<ui> *external = nullptr;
  if (wholeRootBranch) {
    ccrBranchX.assign(firstRow.begin(), firstRow.begin() + forwardOffset);
    external = &ccrBranchX;
  } else {
    ccrCommon.assign(firstRow.begin(), firstRow.end());
    for (size_t i = 1; i < M.size() && !ccrCommon.empty(); ++i) {
      intersectInto(ccrCommonScratch, ccrCommon, adjacentVertices(M[i]));
      ccrCommon.swap(ccrCommonScratch);
    }
  }
  addCcrMetric(ccrFullCalls, 1);
  if (Q.size() <= 4) {
    runSmallCcrBranch(M, Q,
                      wholeRootBranch ? *external : ccrCommon,
                      nullptr);
    return;
  }
  if (!wholeRootBranch)
    setDiffInto(ccrBranchX, ccrCommon, Q);
  prepareCcrBranch(Q, ccrBranchX, true);
  if (ccrWordCount == 1) {
    enumeratePreparedOneWordBranch(M);
    return;
  }
  if (ccrWordCount == 2 && Q.size() >= 68) {
    enumeratePreparedTwoWordBranch(M);
    return;
  }
  ccrBranchR.assign(M.begin(), M.end());
  enumerateAllPureBranchRecursive(ccrBranchR, ccrStates[0], true, 0);
}

void ReorderSib::enumerateAllPureBranchRecursive(
    vector<ui> &R, CcrBitState &state, bool fromP, size_t depth) {
  addCcrMetric(ccrFullStates, 1);
  if (R.size() + state.coreCount + state.residualCount < minCliqueSize)
    return;

  const size_t entryRSize = R.size();
  normalizeCcrState(R, state, fromP);

  if (fromP && R.size() + state.coreCount >= minCliqueSize &&
      ccrCoreUnionIsMaximal(state)) {
    ccrCliqueScratch.assign(R.begin(), R.end());
    ccrForEachBit(state.core, [&](ui local) {
      ccrCliqueScratch.push_back(ccrQVertices[local]);
    });
    recordPureClique(ccrCliqueScratch);
  }
  if (state.residualCount == 0) {
    R.resize(entryRSize);
    return;
  }

  pruneCcrExcludedWithoutResidualNeighbors(state);
  if (state.residualCount == 1) {
    const ui local = ccrFirstBit(state.residual);
    CcrBitState &child = ccrStates[depth + 1];
    buildCcrChildState(child, state, local);
    R.push_back(ccrQVertices[local]);
    enumerateAllPureBranchRecursive(R, child, true, depth + 1);
    R.resize(entryRSize);
    return;
  }

  const size_t pivot = selectCcrPivot(state);
  buildCcrBranchRoots(state, pivot);
  for (ui local : state.branchRoots) {
    const bool selectedFromP = ccrBitIsSet(state.residual, local);
    if (!selectedFromP && !ccrBitIsSet(state.core, local))
      throw logic_error("CCRMCE full branch root was already removed");

    CcrBitState &child = ccrStates[depth + 1];
    buildCcrChildState(child, state, local);

    R.push_back(ccrQVertices[local]);
    enumerateAllPureBranchRecursive(R, child, selectedFromP,
                                    depth + 1);
    R.pop_back();

    if (selectedFromP) {
      ccrClearBit(state.residual, local);
      --state.residualCount;
    } else {
      ccrClearBit(state.core, local);
      --state.coreCount;
    }
    ccrSetBit(state.excludedCandidates, local);
  }
  R.resize(entryRSize);
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

bool ReorderSib::recordPureClique(vector<ui> &C) {
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
      return false;
    }
  }

  if (storedCliqueCount() > numeric_limits<ui>::max())
    throw overflow_error("materialized clique index exceeds uint32_t");

  const ui cliqueId = static_cast<ui>(storedCliqueCount());
  addCliqueCountOrThrow(cliqueCount, 1);
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
  // A posting is queried only from branches that exceed the conservative
  // direct-CCRMCE threshold. Every such branch keeps its original root in M,
  // so indexing vertices whose own forward root is necessarily direct only
  // wastes insertion time and memory without making a cover discoverable.
  for (ui v : C) {
    const size_t forwardCount =
        adjacentVertices(v).size() - firstForwardNeighbor[v];
    if (forwardCount > kSmallQCcrThreshold)
      cliqueIdsByVertex[v].push_back(cliqueId);
  }
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
  cliqueVertices.clear();
  cliqueOffsets.clear();
  cliqueOffsets.push_back(0);
  emittedHashKeys.clear();
  emittedHashHeads.clear();
  emittedHashNext.clear();
  emittedHashSlotsUsed = 0;
  cliqueIdsByVertex.assign(n, {});
  ccrFindOneCalls = 0;
  ccrFullCalls = 0;
  ccrFindOneStates = 0;
  ccrFullStates = 0;
  ccrCoreExtractions = 0;
  ccrCoreVertices = 0;
  ccrResidualVertices = 0;
  solverBudgetFallbacks = 0;
  solverCapacityFallbacks = 0;
  seedSolverCalls = 0;
  maximumSeedConstraints = 0;

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
    worklist.push_back(std::move(nextBranch));
  };

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
    if (1 + forwardCount < minCliqueSize)
      continue;

    // Cache the wider-route productivity decision once per relevant root.
    // Low-degree roots skip the ratio arithmetic, while roots with many
    // sibling branches avoid recomputing the same decision on every pop.
    bool useAdaptiveDirectRoute = false;
    if (forwardCount > kSmallQCcrThreshold) {
      const ull completedRoots = root;
      const ull minimumProductiveCliques =
          (completedRoots +
           kAdaptiveDirectMinCliquesPerRootDenominator - 1) /
          kAdaptiveDirectMinCliquesPerRootDenominator;
      useAdaptiveDirectRoute =
          completedRoots < kAdaptiveDirectWarmupRoots ||
          cliqueCount >= minimumProductiveCliques;
    }
    // Every clique below this branch contains root and only vertices larger
    // than root. Different roots therefore cannot produce equal cliques.
    // Keep global clique storage and cover postings, but restart only the
    // exact duplicate hash at a small root-local capacity.
    constexpr size_t kInitialRootHashCapacity = 16;
    emittedHashKeys.resize(kInitialRootHashCapacity);
    emittedHashHeads.resize(kInitialRootHashCapacity);
    fill(emittedHashHeads.begin(), emittedHashHeads.end(),
         numeric_limits<size_t>::max());
    emittedHashSlotsUsed = 0;

    PureBranch rootBranch = acquireBranch();
    rootBranch.mustin.clear();
    rootBranch.expandTo.clear();
    rootBranch.mustin.push_back(root);
    rootBranch.expandTo.assign(forwardBegin, neighbors.end());
    pushBranch(std::move(rootBranch));
    const auto processRoot = [&](auto directThresholdTag) {
      constexpr size_t directCcrThreshold =
          decltype(directThresholdTag)::value;
      while (!worklist.empty()) {
      PureBranch branch = std::move(worklist.back());
      worklist.pop_back();

      if (branch.mustin.size() + branch.expandTo.size() < minCliqueSize) {
        recycleBranch(std::move(branch));
        continue;
      }

      if (branch.expandTo.size() <= directCcrThreshold) {
        // Small branches use the same exhaustive CCRMCE engine. Its
        // one-word state is cheaper than entering a separate subset kernel
        // and keeps every enumeration route on one implementation.
        enumerateAllPureBranch(branch.mustin, branch.expandTo);
        recycleBranch(std::move(branch));
        continue;
      }

      if (!collectAllCoveringCliques(branch.mustin,
                                     coveringCliqueScratch)) {
        enumerateAllPureBranch(branch.mustin, branch.expandTo);
        recycleBranch(std::move(branch));
        continue;
      }
      const vector<ui> &covers = coveringCliqueScratch;
      if (covers.empty()) {
        vector<ui> found;
        if (findOnePure(branch.mustin, branch.expandTo, found)) {
          recordPureClique(found);
          // The unchanged branch retains every unseen target. On its next pop,
          // the clique just recorded necessarily covers its must-in set.
          pushBranch(std::move(branch));
        } else {
          recycleBranch(std::move(branch));
        }
        continue;
      }

      bool usePivotFallback = false;
      vector<vector<ui>> seeds = generateExactSiblingSets(
          branch.expandTo, covers, &usePivotFallback);
      if (usePivotFallback) {
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
        if (nextBranch.mustin.size() + nextBranch.expandTo.size() <
            minCliqueSize) {
          recycleBranch(std::move(nextBranch));
          continue;
        }
        pushBranch(std::move(nextBranch));
      }
      recycleBranch(std::move(branch));
      }
    };
    if (useAdaptiveDirectRoute)
      processRoot(integral_constant<size_t, kAdaptiveDirectQThreshold>{});
    else
      processRoot(integral_constant<size_t, kSmallQCcrThreshold>{});
  }

  if (static_cast<ull>(storedCliqueCount()) != cliqueCount)
    throw logic_error("stored Pure clique count differs from numeric count");

  cout << "reorder.cliques=" << cliqueCount << '\n'
       << "reorder.stored_cliques=" << storedCliqueCount() << '\n'
       << "reorder.minimum_clique_size=" << minCliqueSize << '\n'
       << "reorder.budget=";
  if (solverWorkBudgetEnabled)
    cout << solverWorkBudget;
  else
    cout << "unlimited";
  cout << '\n'
       << "reorder.ccr.enabled=1\n"
       << "reorder.config.et1=" << kEt1Enabled << '\n'
       << "reorder.config.et2=" << kEt2Enabled << '\n'
       << "reorder.config.et3=" << kEt3Enabled << '\n'
       << "reorder.config.hitset_capacity=";
  if (kHitsetDynamic)
    cout << "dynamic";
  else
    cout << kHitsetCapacity;
  cout << '\n'
       << "reorder.config.hitset_dynamic=" << kHitsetDynamic << '\n'
       << "reorder.config.small_q_ccr_threshold=" << kSmallQCcrThreshold << '\n'
       << "reorder.config.adaptive_direct_q_threshold="
       << kAdaptiveDirectQThreshold << '\n'
       << "reorder.config.pruning.normalization=" << kPruneNormalization << '\n'
       << "reorder.config.pruning.subsumption=" << kPruneSubsumption << '\n'
       << "reorder.config.pruning.unit=" << kPruneUnit << '\n'
       << "reorder.config.pruning.usefulness=" << kPruneUsefulness << '\n'
       << "reorder.config.pruning.antichain=" << kPruneAntichain << '\n'
       << "reorder.config.pruning.fail_first=" << kPruneFailFirst << '\n'
       << "reorder.config.pruning.zero_coverage=" << kPruneZeroCoverage << '\n'
       << "reorder.seed_solver_calls=" << seedSolverCalls << '\n'
       << "reorder.maximum_seed_constraints=" << maximumSeedConstraints << '\n'
       << "reorder.budget_fallbacks=" << solverBudgetFallbacks << '\n'
       << "reorder.capacity_fallbacks=" << solverCapacityFallbacks << '\n'
       << "reorder.ccr.findone_calls=" << ccrFindOneCalls << '\n'
       << "reorder.ccr.full_calls=" << ccrFullCalls << '\n'
       << "reorder.ccr.findone_states=" << ccrFindOneStates << '\n'
       << "reorder.ccr.full_states=" << ccrFullStates << '\n'
       << "reorder.ccr.core_extractions=" << ccrCoreExtractions << '\n'
       << "reorder.ccr.core_vertices=" << ccrCoreVertices << '\n'
       << "reorder.ccr.residual_vertices=" << ccrResidualVertices << '\n';
}
