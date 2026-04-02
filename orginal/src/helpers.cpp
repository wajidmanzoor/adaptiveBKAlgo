#include "../inc/helpers.h"
#include <chrono>
#include <numeric>

AdjMatBK::AdjMatBK(Graph &g) {
  n = g.n;
  cliqueCount = 0;
  adjMatrix.resize(n, vector<bool>(n, false));

  // Initialize adjacency matrix from graph's adjacency list
  for (ui u = 0; u < n; u++) {
    for (ui idx = g.offset[u]; idx < g.offset[u + 1]; idx++) {
      ui v = g.neighbors[idx];
      adjMatrix[u][v] = true;
      adjMatrix[v][u] = true; // Undirected graph
    }
  }
}

vector<bool> AdjMatBK::intersect(const vector<bool> &set1,
                                 const vector<bool> &set2) {
  vector<bool> result(n, false);
  for (ui i = 0; i < n; i++) {
    result[i] = set1[i] && set2[i];
  }
  return result;
}

vector<bool> AdjMatBK::setUnion(const vector<bool> &set1,
                                const vector<bool> &set2) {
  vector<bool> result(n, false);
  for (ui i = 0; i < n; i++) {
    result[i] = set1[i] || set2[i];
  }
  return result;
}

vector<bool> AdjMatBK::setDifference(const vector<bool> &set1,
                                     const vector<bool> &set2) {
  vector<bool> result(n, false);
  for (ui i = 0; i < n; i++) {
    result[i] = set1[i] && !set2[i];
  }
  return result;
}

bool AdjMatBK::isEmpty(const vector<bool> &set) {
  for (bool val : set) {
    if (val)
      return false;
  }
  return true;
}

void AdjMatBK::printSet(const vector<bool> &set, const string &name) {
  cout << name << ": { ";
  for (ui i = 0; i < n; i++) {
    if (set[i]) {
      cout << i << " ";
    }
  }
  cout << "}" << endl;
}

void AdjMatBK::bronKerboschRecursive(vector<bool> &R, vector<bool> &P,
                                     vector<bool> &X) {
  if (isEmpty(P) && isEmpty(X)) {
    // Found a maximal clique
    cliqueCount++;
    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui i = 0; i < n; i++) {
        if (R[i]) {
          cout << i << " ";
        }
      }
      cout << "}" << endl;
    }
    return;
  }

  vector<bool> P_copy = P; // Copy of P to iterate over

  for (ui v = 0; v < n; v++) {
    if (P_copy[v]) {

      vector<bool> new_R = R;
      new_R[v] = true;

      // Create P ∩ N(v) and X ∩ N(v)
      vector<bool> new_P = intersect(P, adjMatrix[v]);
      vector<bool> new_X = intersect(X, adjMatrix[v]);

      // Recursive call
      bronKerboschRecursive(new_R, new_P, new_X);

      // Move v from P to X
      P[v] = false;
      X[v] = true;
    }
  }
}

void AdjMatBK::findAllMaximalCliques() {
  vector<bool> R(n, false); // Current clique
  vector<bool> P(n, true);  // Potential candidates
  vector<bool> X(n, false); // Already processed
  cliqueCount = 0;
  bronKerboschRecursive(R, P, X);
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
}

// Adjacency List based Bron-Kerbosch
AdjListBK::AdjListBK(Graph &g) {
  n = g.n;
  adjList.resize(n);

  // Fill adjacency lists from graph
  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < n) {
        adjList[i].push_back(neighbor);
      }
    }
    // Sort for efficient intersection operations
    sort(adjList[i].begin(), adjList[i].end());
  }
}

vector<ui> AdjListBK::intersect(const vector<ui> &set1,
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

bool AdjListBK::isEmpty(const vector<ui> &set) { return set.empty(); }

void AdjListBK::bronKerboschRecursive(vector<ui> &R, vector<ui> &P,
                                      vector<ui> &X) {
  if (isEmpty(P) && isEmpty(X)) {
    // Found a maximal clique
    cliqueCount++;
    if (R.size() > 2) {
      ofstream outfile("bk_adj_maximal_cliques.txt", std::ios::app);
      if (outfile.is_open()) {
        outfile << "Maximal Clique Found: { ";
        for (ui x : R)
          outfile << x << " ";
        outfile << "}\n";
      }
    }
    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui v : R) {
        cout << v << " ";
      }
      cout << "}" << endl;
    }
    return;
  }

  vector<ui> P_copy = P; // Copy of P to iterate over

  for (ui v : P_copy) {
    vector<ui> new_R = R;
    new_R.push_back(v);

    // Create P ∩ N(v) and X ∩ N(v)
    vector<ui> new_P = intersect(P, adjList[v]);
    vector<ui> new_X = intersect(X, adjList[v]);

    // Recursive call
    bronKerboschRecursive(new_R, new_P, new_X);

    // Move v from P to X
    P.erase(find(P.begin(), P.end(), v));
    X.push_back(v);
  }
}

void AdjListBK::findAllMaximalCliques() {
  // Initialize sets
  vector<ui> R; // Current clique (empty)
  vector<ui> P; // All vertices as candidates
  vector<ui> X; // Excluded set (empty)

  // Fill P with all vertices
  for (ui i = 0; i < n; i++) {
    P.push_back(i);
  }

  cliqueCount = 0;
  bronKerboschRecursive(R, P, X);

  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
}

// PivotBK Implementation
PivotBK::PivotBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  // Fill adjacency lists from graph
  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < n) {
        adjList[i].push_back(neighbor);
      }
    }
    // Sort for efficient intersection operations
    sort(adjList[i].begin(), adjList[i].end());
  }
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
  if (isEmpty(P) && isEmpty(X)) {
    if (R.size() <= 2)
      return;

    // Found a maximal clique
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());

    ofstream outfile("bk_pivot_maximal_cliques.txt", std::ios::app);
    if (outfile.is_open()) {
      outfile << "Maximal Clique Found: { ";
      for (ui x : R)
        outfile << x << " ";
      outfile << "}\n";
    }
    foundCliques.push_back(R);

    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui v : R) {
        cout << v << " ";
      }
      cout << "}" << endl;
    }
    return;
  }

  // P empty but X non-empty: R is NOT maximal, prune
  if (isEmpty(P))
    return;

  // Choose pivot from P ∪ X
  ui pivot = choosePivot(P, X);

#if debug
  if (isPSubsetOfFoundClique(P)) {
    redendantChecks.push_back(P);
    redundancy++;
  }
#endif

  // Snapshot P \ N(pivot) before modifying P
  vector<ui> candidates;
  candidates.reserve(P.size());
  for (ui v : P) {
    if (!isConnected(v, pivot))
      candidates.push_back(v);
  }

  for (ui v : candidates) {
    R.push_back(v);
    vector<ui> new_P = intersect(P, adjList[v]);
    vector<ui> new_X = intersect(X, adjList[v]);

    bronKerboschRecursive(R, new_P, new_X);

    R.pop_back(); // O(1) backtrack — no copy needed

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
  foundCliques.clear();

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

static void printVec(const string &label, const vector<ui> &v);
static void printForest(const string &label, const vector<vector<ui>> &mustin,
                        const vector<vector<ui>> &expandTo);

// ── print helpers
// ─────────────────────────────────────────────────────────────
static void printVec(const string &label, const vector<ui> &v) {
  cout << label << "{ ";
  for (ui x : v)
    cout << x << " ";
  cout << "}" << endl;
}

static void printForest(const string &label, const vector<vector<ui>> &mustin,
                        const vector<vector<ui>> &expandTo) {
  cout << label << endl;
  for (ui i = 0; i < mustin.size(); i++) {
    cout << "    [" << i << "]  mustin=[ ";
    for (ui v : mustin[i])
      cout << v << " ";
    cout << "]  expandTo=[ ";
    for (ui v : expandTo[i])
      cout << v << " ";
    cout << "]" << endl;
  }
}

// ═════════════════════════════════════════════════════════════════════════════
// ReorderNew
// ═════════════════════════════════════════════════════════════════════════════

ReorderNew::ReorderNew(Graph &g) {
  n = g.n;
  adjList.resize(n);
  adjList2.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  cliquesByVertex.resize(n);

  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui nb = g.neighbors[j];
      adjList[i].push_back(nb);
      if (nb > i)
        adjList2[i].push_back(nb);
    }
    sort(adjList[i].begin(), adjList[i].end());
    sort(adjList2[i].begin(), adjList2[i].end());
  }
}

// ─── set helpers ─────────────────────────────────────────────────────────────

vector<ui> ReorderNew::intersect(vector<ui> A, vector<ui> B) {
  vector<char> seen(n, 0);
  for (ui x : A)
    seen[x] = 1;
  vector<ui> C;
  for (ui y : B)
    if (seen[y])
      C.push_back(y);
  return C;
}

vector<ui> ReorderNew::setDiff(vector<ui> A, vector<ui> B) {
  vector<char> seen(n, 0);
  for (ui x : B)
    seen[x] = 1;
  vector<ui> C;
  for (ui y : A)
    if (!seen[y])
      C.push_back(y);
  return C;
}

vector<ui> ReorderNew::unionSet(vector<ui> A, vector<ui> B) {
  vector<char> seen(n, 0);
  vector<ui> U;
  for (ui x : A) {
    if (!seen[x]) {
      seen[x] = 1;
      U.push_back(x);
    }
  }
  for (ui y : B) {
    if (!seen[y]) {
      seen[y] = 1;
      U.push_back(y);
    }
  }
  return U;
}

// Intersection of adjList[v] for all v in mustin — vertices adjacent to all of mustin
vector<ui> ReorderNew::commonNeighbors(const vector<ui> &mustin) {
  if (mustin.empty())
    return {};
  vector<ui> result = adjList[mustin[0]];
  for (ui i = 1; i < (ui)mustin.size(); i++)
    result = intersect(result, adjList[mustin[i]]);
  return result;
}

// ─── applyEffect ─────────────────────────────────────────────────────────────
// Apply the effect of found clique C on tree (mustinJ, expandToJ).
//
// Rule 1 — Skip:       mustinJ ⊆ C  →  return true
// Restriction has two cases depending on whether this tree's root == min(C):
//
//  min(M) == min(C):  remove ALL C\M from expandTo (this tree found C; its
//                     remaining sub-trees would re-find C or non-maximal subsets).
//                     Skip if expandTo becomes empty.
//
//  min(M) >  min(C):  remove only {v ∈ C : v < min(M)} from expandTo.
//                     These are lower-index vertices that were added by
//                     enrichment; keeping them would reconstruct C via a
//                     reverse path.  Higher-index C members stay so paths to
//                     genuinely new cliques are not blocked.
//
// Rule 3 — Enrich: mustinJ ∩ C = ∅ → expandToJ += (C ∩ commonNeighbors(mustinJ))

bool ReorderNew::applyEffect(vector<ui> &mustinJ, vector<ui> &expandToJ,
                              const vector<ui> &C) {
  vector<char> inC(n, 0);
  for (ui v : C)
    inC[v] = 1;

  // Check mustinJ ⊆ C
  bool allInC = true;
  for (ui v : mustinJ)
    if (!inC[v]) {
      allInC = false;
      break;
    }

  if (allInC) {
    ui minM = *min_element(mustinJ.begin(), mustinJ.end());
    ui minC = *min_element(C.begin(), C.end());

    vector<char> inMustin(n, 0);
    for (ui v : mustinJ)
      inMustin[v] = 1;

    if (minM == minC) {
      // Full restriction: remove all C\M — would re-find C from this root
      vector<ui> toRemove;
      for (ui v : C)
        if (!inMustin[v])
          toRemove.push_back(v);
      expandToJ = setDiff(expandToJ, toRemove);
      return expandToJ.empty();
    } else {
      // Partial restriction: remove enriched lower-index C members
      vector<ui> toRemove;
      for (ui v : C)
        if (v < minM && !inMustin[v])
          toRemove.push_back(v);
      if (!toRemove.empty())
        expandToJ = setDiff(expandToJ, toRemove);
      // If every remaining candidate is still inside C, all DFS paths lead to
      // subsets of C (non-maximal) — skip the tree entirely
      for (ui v : expandToJ)
        if (!inC[v])
          return false; // at least one non-C candidate exists
      return expandToJ.empty() ? false : true; // non-empty but all-in-C → skip
    }
  }

  bool hasOverlap = false;
  for (ui v : mustinJ)
    if (inC[v]) {
      hasOverlap = true;
      break;
    }

  if (!hasOverlap) {
    // Enrichment: add C members adjacent to all of mustinJ
    vector<ui> cn = commonNeighbors(mustinJ);
    vector<ui> toAdd = intersect(vector<ui>(C.begin(), C.end()), cn);
    expandToJ = unionSet(expandToJ, toAdd);
  }
  // Mixed case: no action
  return false;
}

// ─── retroRestrict ───────────────────────────────────────────────────────────
// Apply restriction-only from ALL cliques found so far that touch mustinI.
// Handles cliques found by earlier sibling rCalls after this tree's expandTo
// was last set. Returns true if the tree should be skipped.

bool ReorderNew::retroRestrict(vector<ui> &mustinI, vector<ui> &expandToI) {
  for (ui v : mustinI) {
    for (ui cliqueId : cliquesByVertex[v]) {
      const vector<ui> &C = allCliques[cliqueId];
      vector<char> inC(n, 0);
      for (ui x : C)
        inC[x] = 1;

      // Only restrict when mustinI ⊆ C
      bool allIn = true;
      for (ui x : mustinI)
        if (!inC[x]) {
          allIn = false;
          break;
        }
      if (!allIn)
        continue;

      ui minM = *min_element(mustinI.begin(), mustinI.end());
      ui minC = *min_element(C.begin(), C.end());

      vector<char> inMustin(n, 0);
      for (ui x : mustinI)
        inMustin[x] = 1;

      if (minM == minC) {
        // Full restriction
        vector<ui> toRemove;
        for (ui x : C)
          if (!inMustin[x])
            toRemove.push_back(x);
        expandToI = setDiff(expandToI, toRemove);
        if (expandToI.empty())
          return true;
      } else {
        // Partial restriction: only lower-index enriched C members
        vector<ui> toRemove;
        for (ui x : C)
          if (x < minM && !inMustin[x])
            toRemove.push_back(x);
        if (!toRemove.empty())
          expandToI = setDiff(expandToI, toRemove);
        // If every remaining candidate is in C, all paths lead to subsets of C
        bool allInC = true;
        for (ui x : expandToI)
          if (!inC[x]) {
            allInC = false;
            break;
          }
        if (allInC && !expandToI.empty())
          return true;
      }
    }
  }
  return false;
}

// ─── rCall ───────────────────────────────────────────────────────────────────

void ReorderNew::rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo,
                       ui level) {
  for (ui i = 0; i < (ui)mustin.size(); i++) {
    if (retroRestrict(mustin[i], expandTo[i]))
      continue;

    vector<ui> R = mustin[i];
    vector<ui> Q = expandTo[i];
    bool done = false;
    enumerate(R, Q, mustin, expandTo, i, level, done);
    if (done)
      break;
  }
}

// ─── enumerate ───────────────────────────────────────────────────────────────

void ReorderNew::enumerate(vector<ui> &R, vector<ui> &Q,
                            vector<vector<ui>> &mustin,
                            vector<vector<ui>> &expandTo, ui treeIndex,
                            ui level, bool &done) {
  checksCount++;

  // ── Base case: no candidates left ────────────────────────────────────────
  if (Q.empty()) {
    if ((ui)R.size() > 2) {
      // Record clique
      cliqueCount++;
      if ((ui)R.size() > maxCliqueSize)
        maxCliqueSize = (ui)R.size();
      ui cliqueIdx = (ui)allCliques.size();
      allCliques.push_back(R);
      for (ui v : R)
        cliquesByVertex[v].push_back(cliqueIdx);

      // Apply effect of R on trees treeIndex..end, then spawn sub-rCalls
      vector<bool> skip(mustin.size(), false);
      for (ui j = treeIndex; j < (ui)mustin.size(); j++)
        skip[j] = applyEffect(mustin[j], expandTo[j], R);

      for (ui j = treeIndex; j < (ui)mustin.size(); j++) {
        if (skip[j])
          continue;

        // Build sub-forest for tree j
        vector<vector<ui>> subMustin;
        vector<vector<ui>> subExpandTo;

        if (expandTo[j].empty()) {
          // No expansion: mustin[j] itself is a maximal clique candidate
          if ((ui)mustin[j].size() > 2) {
            subMustin.push_back(mustin[j]);
            subExpandTo.push_back({});
          }
        } else {
          for (ui v : expandTo[j]) {
            vector<ui> newMustin = mustin[j];
            newMustin.push_back(v);
            vector<ui> newET = intersect(expandTo[j], adjList2[v]);
            subMustin.push_back(newMustin);
            subExpandTo.push_back(newET);
          }
        }

        if (!subMustin.empty())
          rCall(subMustin, subExpandTo, level + 1);
      }
      done = true;
    }
    return;
  }

  // ── Recursive case: expand R by each candidate in Q ──────────────────────
  for (ui v : Q) {
    R.push_back(v);
    vector<ui> Qp = intersect(Q, adjList2[v]);
    enumerate(R, Qp, mustin, expandTo, treeIndex, level, done);
    R.pop_back();
    if (done)
      return;
  }
}

// ─── findAllMaximalCliques ────────────────────────────────────────────────────

void ReorderNew::findAllMaximalCliques() {
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  allCliques.clear();
  for (ui v = 0; v < n; v++)
    cliquesByVertex[v].clear();

  vector<vector<ui>> mustin;
  vector<vector<ui>> expandTo;
  for (ui v = 0; v < n; v++) {
    mustin.push_back({v});
    expandTo.push_back(adjList2[v]);
  }

  auto t0 = chrono::high_resolution_clock::now();
  rCall(mustin, expandTo, 0);
  auto t1 = chrono::high_resolution_clock::now();
  double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << "ReorderNew: cliques=" << cliqueCount
       << "  maxSize=" << maxCliqueSize
       << "  checks=" << checksCount
       << "  time=" << ms << " ms" << endl;
}

// ─────────────────────────────────────────────────────────────────────────────
Reorder::Reorder(Graph &g) {
  graph = g;
  n = g.n;
  adjList.resize(n);
  adjList2.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  // status array removed
  cliquesByVertex.resize(n);

  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < n) {
        adjList[i].push_back(neighbor);
      }
      if ((neighbor < n) && (neighbor > i)) {
        adjList2[i].push_back(neighbor);
      }
    }
    sort(adjList[i].begin(), adjList[i].end());
    sort(adjList2[i].begin(), adjList2[i].end());
  }
}

void Reorder::findAllMaximalCliques() {

  cliqueCount = 0;
  maxCliqueSize = 0;
  allCliques.clear();
  for (ui v = 0; v < n; v++)
    cliquesByVertex[v].clear();

  vector<vector<ui>> mustin;
  vector<vector<ui>> expandTo;

  for (ui v = 0; v < n; v++) {
    mustin.push_back({v});
    expandTo.push_back(adjList2[v]);
  }

  cout << "\n[INIT] Initial forest:" << endl;
  printForest("", mustin, expandTo);

  rCall(mustin, expandTo, 0, 0);

  cout << "\n╔══════════════════════════════════════════════╗" << endl;
  cout << "║  DONE  cliques=" << cliqueCount << "  maxSize=" << maxCliqueSize
       << endl;
  cout << "╚══════════════════════════════════════════════╝" << endl;

  cout << "\n[RESULT] allCliques (" << allCliques.size() << " total):" << endl;
  for (ui i = 0; i < allCliques.size(); i++) {
    cout << "  [" << i << "] { ";
    for (ui v : allCliques[i])
      cout << v << " ";
    cout << "}" << endl;
  }

  cout << "\n[RESULT] cliquesByVertex:" << endl;
  for (ui v = 0; v < n; v++) {
    if (!cliquesByVertex[v].empty()) {
      cout << "  vertex " << v << " -> clique indices { ";
      for (ui idx : cliquesByVertex[v])
        cout << idx << " ";
      cout << "}" << endl;
    }
  }
}

vector<ui> Reorder::compliment(vector<ui> &vector1) {
  vector<char> seen(n, 0);
  for (ui x : vector1)
    seen[x] = 1;
  vector<ui> C;
  for (ui i = 0; i < n; i++)
    if (!seen[i])
      C.push_back(i);
  return C;
}

// ─────────────────────────────────────────────────────────────────────────────
void Reorder::rCall(vector<vector<ui>> &mustin, vector<vector<ui>> &expandTo,
                    ui level, ui enlevel) {

  string indent(level * 2, ' ');

  cout << "\n" << indent << "┌─[rCall level=" << level << "]" << endl;
  printForest(indent + "│  Forest in:", mustin, expandTo);

  if (mustin.empty()) {
    cout << indent << "│  (empty forest — return)" << endl;
    cout << indent << "└─[rCall level=" << level << " done]" << endl;
    return;
  }

  vector<ui> Q, R;

  bool moveToNext = true;
  ui index = 0;

  for (ui i = 0; i < mustin.size(); i++) {

    if (!moveToNext) {
      cout << indent << "│  moveToNext=false — stopping rCall loop" << endl;
      break;
    }

    if (mustin[i].empty()) {
      cout << indent << "│  [tree " << i << "] empty mustin — skip" << endl;
      continue;
    }

    // Sibling effect:

    cout << " Sibling effect for tree[" << i << "]" << endl;

    ui skip = 0;
    for (ui j = 0; j < mustin[i].size() - 1; j++) {
      ui v = mustin[i][j];
      for (ui cliqueId : cliquesByVertex[v]) {
        cout << "- Vertex " << v << " affects clique " << cliqueId << endl;
        vector<ui> clique = allCliques[cliqueId];
        cout << "clique ";
        for (ui x : clique)
          cout << x << "";
        cout << endl;
        vector<ui> comp = compliment(clique);
        if (find(comp.begin(), comp.end(), mustin[i].back()) != comp.end()) {
          skip = 1;
          break;
        }
      }
    }
    if (skip) {
      continue;
    }

    cout << " New ExpandTo ";
    for (ui x : expandTo[i]) {
      cout << x << " ";
    }
    cout << endl;

    R.clear();
    for (ui v : mustin[i])
      R.push_back(v);
    Q = expandTo[i];

    cout << indent << "│" << endl;
    cout << indent << "│  [tree " << i << "] starting enumeration" << endl;
    printVec(indent + "│    R=", R);
    printVec(indent + "│    Q=", Q);

    bool flag = false;
    enemurate(R, Q, mustin, expandTo, moveToNext, flag, index, level, enlevel);
    index++;
  }

  cout << indent << "└─[rCall level=" << level << " done]" << endl;
}

// ─────────────────────────────────────────────────────────────────────────────
void Reorder::enemurate(vector<ui> &R, vector<ui> &Q,
                        vector<vector<ui>> &mustin,
                        vector<vector<ui>> &expandTo, bool &moveToNext,
                        bool &flag, ui index, ui level, ui enlevel) {

  string indent(level * 2 + 4, ' ');

  // ── BASE CASE: Q is empty → maximal clique candidate ─────────────────────
  if (Q.empty()) {

    if (R.size() > 2) {

      // ── CLIQUE FOUND ──────────────────────────────────────────────────────
      cliqueCount++;
      if (R.size() > maxCliqueSize)
        maxCliqueSize = (ui)R.size();

      ui cliqueIdx = (ui)allCliques.size();
      allCliques.push_back(R);
      for (ui v : R)
        cliquesByVertex[v].push_back(cliqueIdx);

      flag = true;
      moveToNext = false;

      cout << indent << "★ CLIQUE #" << cliqueCount << " (idx=" << cliqueIdx
           << "): { ";
      for (ui v : R)
        cout << v << " ";
      cout << "}" << endl;

      vector<ui> skip(mustin.size(), 0);
      for (ui i = index; i < mustin.size(); i++) {
        vector<ui> tree = mustin[i];

        ui len = (ui)tree.size();
        for (ui v : tree)
          if (find(R.begin(), R.end(), v) != R.end())
            len--;

        if (len == 0) {
          skip[i] = 1;
          continue;
        }

        vector<ui> temp2 = unionSet(R, expandTo[i]);

        for (ui v : tree) {
          temp2 = intersect(temp2, adjList[v]);
        }

        expandTo[i] = temp2;
      }

      printForest(indent + "  ", mustin, expandTo);

      for (ui i = 0; i < mustin.size(); i++) {
        if (skip[i]) {
          continue;
        }
        vector<vector<ui>> mustin_;
        vector<vector<ui>> expand_;

        if (expandTo[i].empty()) {
          mustin_.push_back(mustin[i]);
          expand_.push_back(expandTo[i]);
        }

        for (ui v : expandTo[i]) {
          vector<ui> temp = mustin[i];
          temp.push_back(v);
          vector<ui> newET = intersect(expandTo[i], adjList2[v]);

          mustin_.push_back(temp);
          expand_.push_back(newET);
        }

        if (mustin_.empty()) {
          continue;
        }
        rCall(mustin_, expand_, level + 1, 0);
      }

      return;

    } else {

      return;
    }
  }

  // ── RECURSIVE CASE: expand R by each candidate in Q ──────────────────────
  cout << indent << "enumerate(level=" << level << " enlevel=" << enlevel
       << endl;

  for (ui v : Q) {

    if (flag) {
      return;
    }

    R.push_back(v);
    vector<ui> Q_ = intersect(Q, adjList2[v]);

    enemurate(R, Q_, mustin, expandTo, moveToNext, flag, index, level,
              enlevel + 1);

    R.pop_back();
  }
}

vector<ui> Reorder::intersect(vector<ui> vector1, vector<ui> vector2) {
  vector<char> seen(n, 0);
  for (ui x : vector1)
    seen[x] = 1;
  vector<ui> C;
  for (ui y : vector2)
    if (seen[y])
      C.push_back(y);
  return C;
}

vector<ui> Reorder::unionSet(vector<ui> vector1, vector<ui> vector2) {
  vector<char> seen(n, 0);
  vector<ui> U;
  for (ui x : vector1) {
    if (!seen[x]) {
      seen[x] = 1;
      U.push_back(x);
    }
  }
  for (ui y : vector2) {
    if (!seen[y]) {
      seen[y] = 1;
      U.push_back(y);
    }
  }
  return U;
}

vector<ui> Reorder::setDifference(vector<ui> A, vector<ui> B) {
  vector<char> seen(n, 0);
  for (ui x : B)
    seen[x] = 1;
  vector<ui> C;
  for (ui y : A)
    if (!seen[y])
      C.push_back(y);
  return C;
}
