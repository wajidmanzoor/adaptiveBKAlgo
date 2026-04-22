#include "../inc/helpers.h"
#include <chrono>
#include <functional>
#include <numeric>

// ── Degeneracy ordering helpers
// ─────────────────────────────────────────────── Returns peelSeq where
// peelSeq[0] = highest-core vertex, peelSeq[n-1] = lowest. Operates on a copy
// of g.degree so the graph is not modified.
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

// Build adjList (all neighbors) and adjList2 (higher-index neighbors only)
// in the relabeled space defined by perm[old_v] = new_v.
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
PivotBK::PivotBK(Graph &g, DegOrder order) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;

  vector<ui> perm(n);
  if (order == DegOrder::ORIGINAL) {
    for (ui i = 0; i < n; i++)
      perm[i] = i;
  } else {
    vector<ui> peelSeq = computePeelSeq(g);
    if (order == DegOrder::ASCENDING) {
      for (ui i = 0; i < n; i++)
        perm[peelSeq[n - 1 - i]] = i; // low-core → low index
    } else {
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

    // Found a maximal clique — only update counters inside timed section
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());
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

/// Compliment reorder

Reorder::Reorder(Graph &g, DegOrder order) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  cliquesByVertex.resize(n);

  vector<ui> perm(n);
  if (order == DegOrder::ORIGINAL) {
    // Canonical Order.
    for (ui i = 0; i < n; i++)
      perm[i] = i;
  } else {
    vector<ui> peelSeq = computePeelSeq(g);
    if (order == DegOrder::ASCENDING) {

      // ascending core value.
      for (ui i = 0; i < n; i++)
        perm[peelSeq[n - 1 - i]] = i;
    } else { // DESCENDING

      // descending core value.

      for (ui i = 0; i < n; i++)
        perm[peelSeq[i]] = i;
    }
  }

  buildAdjLists(g, perm, adjList, adjList2);
}

vector<ui> Reorder::intersect(vector<ui> A, vector<ui> B) {
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

vector<ui> Reorder::setDiff(vector<ui> A, vector<ui> B) {
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

vector<ui> Reorder::compliment(const vector<ui> &vector1) {
  vector<char> seen(n, 0);
  for (ui x : vector1)
    seen[x] = 1;
  vector<ui> C;
  for (ui i = 0; i < n; i++)
    if (!seen[i])
      C.push_back(i);
  return C;
}

vector<ui> Reorder::unionSet(vector<ui> A, vector<ui> B) {
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
  while (i < A.size()) {
    U.push_back(A[i]);
    i++;
  }
  while (j < B.size()) {
    U.push_back(B[j]);
    j++;
  }
  return U;
}

void Reorder::rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo,
                    ui level) {
  if (debug) {
    for (ui i = 0; i < level; i++) {
      cout << "     ";
      ;
    }
    cout << "Level " << level << ": mustin and expandTo sets:" << endl;
    for (ui i = 0; i < mustin.size(); i++) {
      for (ui i = 0; i < level; i++) {
        cout << "     ";
      }
      cout << "Vertex " << i << ": mustin={ ";
      for (ui v : mustin[i])
        cout << v << " ";
      cout << "}  expandTo={ ";
      for (ui v : expandTo[i])
        cout << v << " ";
      cout << "}" << endl;
    }
  }

  if (level != 0) {
    if (!expandTo[0].empty()) {

      vector<vector<ui>> tempmustin = mustin;
      vector<vector<ui>> tempexpandTo = expandTo;

      vector<ui> temp;
      bool oneFullFound = false;

      for (ui v : mustin[0]) {
        for (ui cId : cliquesByVertex[v]) {
          if (foundLevel[cId] < level - 1)
            continue;
          bool containsAllMustin = true;
          for (ui mv : tempmustin[0]) {
            if (!binary_search(allCliques[cId].begin(), allCliques[cId].end(),
                               mv)) {
              containsAllMustin = false;
              break;
            }
          }
          if (containsAllMustin) {
            oneFullFound = true;
            temp = unionSet(
                temp, intersect(tempexpandTo[0], compliment(allCliques[cId])));
          }
        }
      }

      if (!oneFullFound)
        temp = tempexpandTo[0];

      mustin.clear();
      expandTo.clear();
      for (ui v : temp) {
        vector<ui> newMustin = tempmustin[0];
        newMustin.push_back(v);
        mustin.push_back(newMustin);
        expandTo.push_back(intersect(tempexpandTo[0], adjList[v]));
      }
      if (debug) {
        for (ui i = 0; i < level; i++) {
          cout << "     ";
        }
        cout << "Level " << level << ": After Sibling Effect:" << endl;
        for (ui i = 0; i < mustin.size(); i++) {
          for (ui j = 0; j < level; j++) {
            cout << "     ";
          }
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
  }
  for (ui i = 0; i < (ui)mustin.size(); i++) {
    vector<ui> R = mustin[i];
    vector<ui> Q = expandTo[i];
    bool done = false;
    enumerate(R, Q, mustin, expandTo, i, level, done);
    if (done)
      break;
  }
}

void Reorder::enumerate(vector<ui> &R, vector<ui> &Q,
                        vector<vector<ui>> &mustin,
                        vector<vector<ui>> &expandTo, ui treeIndex, ui level,
                        bool &done) {
  checksCount++;

  if (debug) {
    for (ui i = 0; i < level; i++) {
      cout << "     ";
      ;
    }
    cout << "Level " << level << ": Checking R={ ";
    for (ui v : R)
      cout << v << " ";
    cout << "}  Q={ ";
    for (ui v : Q)
      cout << v << " ";
    cout << "}" << endl;
  }

  if (Q.empty()) {
    if ((ui)R.size() > 2) {
      vector<ui> C = R;
      sort(C.begin(), C.end());
      cliqueCount++;
      if (debug) {
        for (ui i = 0; i < level; i++) {
          cout << "   ";
        }
        cout << "Maximal Clique Found: { ";
        for (ui v : C)
          cout << v << " ";
        cout << "}" << endl;
      }
      if ((ui)C.size() > maxCliqueSize)
        maxCliqueSize = (ui)C.size();
      ui cliqueIdx = (ui)allCliques.size();
      allCliques.push_back(C);
      foundLevel.push_back(level);
      for (ui v : C)
        cliquesByVertex[v].push_back(cliqueIdx);

      vector<ui> comp = compliment(C);
      done = true;

      if (level == 0) {
        cout << "here" << endl;
        vector<vector<ui>> newMustin;
        vector<vector<ui>> newExpandTo;
        for (ui i = treeIndex; i < (ui)mustin.size(); i++) {
          if (find(C.begin(), C.end(), mustin[i].back()) != C.end())
            continue;
          newMustin.push_back(mustin[i]);
          newExpandTo.push_back(
              intersect(setDiff(adjList[mustin[i].back()], mustin[i]),
                        unionSet(C, expandTo[i])));
          if (debug) {
            for (ui i = 0; i < level; i++) {
              cout << "   ";
            }
            cout << "Level " << level << ": After Reorder:" << endl;
            for (ui i = 0; i < newMustin.size(); i++) {
              for (ui j = 0; j < level; j++) {
                cout << "     ";
              }
              cout << "Vertex " << i << ": mustin={ ";
              for (ui v : newMustin[i])
                cout << v << " ";
              cout << "}  expandTo={ ";
              for (ui v : newExpandTo[i])
                cout << v << " ";
              cout << "}" << endl;
            }
          }
          rCall(std::move(newMustin), std::move(newExpandTo), level + 1);
        }

        return;

      } else {
        // Same as level 0: generate seeds for all complement vertices adjacent
        // to C
        vector<vector<ui>> newMustin;
        vector<vector<ui>> newExpandTo;
        cout << "new here" << endl;

        for (ui i = treeIndex; i < (ui)mustin.size(); i++) {
          /*cout << "Tree Index " << i << " mustin[i].back() " <<
             mustin[i].back()
               << endl;*/
          if (find(C.begin(), C.end(), mustin[i].back()) != C.end())
            continue;
          vector<ui> newExp =
              intersect(setDiff(adjList[mustin[i].back()], mustin[i]),
                        unionSet(C, expandTo[i]));
          newMustin.push_back(mustin[i]);
          newExpandTo.push_back(newExp);
          if (debug) {
            for (ui i = 0; i < level; i++) {
              cout << "     ";
            }
            cout << "Level " << level << ": After Reorder:" << endl;
            for (ui i = 0; i < newMustin.size(); i++) {
              for (ui j = 0; j < level; j++) {
                cout << "     ";
              }
              cout << "Vertex " << i << ": mustin={ ";
              for (ui v : newMustin[i])
                cout << v << " ";
              cout << "}  expandTo={ ";
              for (ui v : newExpandTo[i])
                cout << v << " ";
              cout << "}" << endl;
            }
          }
          rCall(std::move(newMustin), std::move(newExpandTo), level + 1);
        }

        return;
      }
    }
  }

  for (ui v : Q) {
    R.push_back(v);
    vector<ui> Qp = intersect(Q, adjList2[v]);
    enumerate(R, Qp, mustin, expandTo, treeIndex, level, done);
    R.pop_back();
    if (done)
      return;
  }
}

void Reorder::findAllMaximalCliques() {
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
  rCall(std::move(mustin), std::move(expandTo), 0);
  auto t1 = chrono::high_resolution_clock::now();
  double ms = chrono::duration<double, milli>(t1 - t0).count();

  cout << "ReorderNew: cliques=" << cliqueCount << "  maxSize=" << maxCliqueSize
       << "  checks=" << checksCount << "  time=" << ms << " ms" << endl;
}

// ReorderSib Implementation
ReorderSib::ReorderSib(Graph &g, DegOrder order, SibMethod method) {
  n = g.n;
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  this->method = method;
  cliquesByVertex.resize(n);

  vector<ui> perm(n);
  if (order == DegOrder::ORIGINAL) {
    for (ui i = 0; i < n; i++)
      perm[i] = i;
  } else {
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

vector<ui> ReorderSib::intersect(vector<ui> A, vector<ui> B) {
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

vector<ui> ReorderSib::setDiff(vector<ui> A, vector<ui> B) {
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

vector<ui> ReorderSib::unionSet(vector<ui> A, vector<ui> B) {
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

bool ReorderSib::isCliqueSet(const vector<ui> &S) {
  for (ui i = 0; i < S.size(); i++) {
    for (ui j = i + 1; j < S.size(); j++) {
      if (!binary_search(adjList[S[i]].begin(), adjList[S[i]].end(), S[j]))
        return false;
    }
  }
  return true;
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

vector<ui> ReorderSib::commonExpand(const vector<ui> &E,
                                    const vector<ui> &S) {
  vector<ui> result = setDiff(E, S);
  for (ui v : S)
    result = intersect(result, adjList[v]);
  return result;
}

vector<ui> ReorderSib::collectCoveringCliques(const vector<ui> &M, ui level) {
  vector<ui> result;
  vector<char> seen(allCliques.size(), 0);
  if (M.empty())
    return result;

  for (ui v : M) {
    for (ui cId : cliquesByVertex[v]) {
      if (cId >= allCliques.size() || seen[cId])
        continue;
      if (foundLevel[cId] < level - 1)
        continue;

      bool containsAllMustin = true;
      for (ui mv : M) {
        if (!binary_search(allCliques[cId].begin(), allCliques[cId].end(),
                           mv)) {
          containsAllMustin = false;
          break;
        }
      }
      if (containsAllMustin) {
        seen[cId] = 1;
        result.push_back(cId);
      }
    }
  }
  return result;
}

vector<vector<ui>> ReorderSib::buildHitSets(const vector<ui> &E,
                                            const vector<ui> &cliqueIds) {
  vector<vector<ui>> hitSets;
  hitSets.reserve(cliqueIds.size());
  for (ui cId : cliqueIds)
    hitSets.push_back(intersect(E, compliment(allCliques[cId])));
  return hitSets;
}

vector<vector<ui>> ReorderSib::singletonBranches(const vector<ui> &E) {
  vector<vector<ui>> branches;
  branches.reserve(E.size());
  for (ui v : E)
    branches.push_back({v});
  return branches;
}

vector<vector<ui>> ReorderSib::minimalBySize(vector<vector<ui>> solutions) {
  for (vector<ui> &S : solutions)
    sort(S.begin(), S.end());
  sort(solutions.begin(), solutions.end());
  solutions.erase(unique(solutions.begin(), solutions.end()), solutions.end());

  if (solutions.empty())
    return solutions;

  size_t bestSize = solutions[0].size();
  for (const vector<ui> &S : solutions)
    bestSize = min(bestSize, S.size());

  vector<vector<ui>> minimal;
  for (const vector<ui> &S : solutions) {
    if (S.size() == bestSize)
      minimal.push_back(S);
  }
  return minimal;
}

vector<vector<ui>>
ReorderSib::minimalByInclusion(vector<vector<ui>> solutions) {
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

vector<vector<ui>> ReorderSib::generateSiblingSets(const vector<ui> &M,
                                                   const vector<ui> &E,
                                                   ui level) {
  vector<ui> cliqueIds = collectCoveringCliques(M, level);
  if (cliqueIds.empty())
    return singletonBranches(E);

  vector<vector<ui>> hitSets = buildHitSets(E, cliqueIds);
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
  }

  return bruteForceBySize(E, hitSets);
}

vector<vector<ui>>
ReorderSib::bruteForceBySize(const vector<ui> &E,
                             const vector<vector<ui>> &hitSets) {
  vector<vector<ui>> solutions;
  vector<ui> current;

  function<void(ui, ui, ui)> choose = [&](ui start, ui remaining,
                                          ui targetSize) {
    if (remaining == 0) {
      if (isCliqueSet(current) && hitsAll(current, hitSets))
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
      choose(i + 1, remaining - 1, targetSize);
      current.pop_back();
    }
  };

  for (ui targetSize = 1; targetSize <= E.size(); targetSize++) {
    choose(0, targetSize, targetSize);
  }
  return minimalByInclusion(solutions);
}

vector<vector<ui>>
ReorderSib::backtrackingBranchBound(const vector<ui> &E,
                                    const vector<vector<ui>> &hitSets) {
  vector<vector<ui>> solutions;
  vector<ui> current;

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

void ReorderSib::backtrackSearch(const vector<ui> &E,
                                 const vector<vector<ui>> &hitSets,
                                 vector<ui> &current, vector<char> &covered,
                                 ui coveredCount, ui start, ui &bestSize,
                                 vector<vector<ui>> &solutions) {
  if (coveredCount == hitSets.size()) {
    if (current.size() < bestSize) {
      bestSize = (ui)current.size();
      solutions.clear();
    }
    if (current.size() == bestSize)
      solutions.push_back(current);
    return;
  }
  if (current.size() >= bestSize)
    return;

  // Pick the first uncovered set, then branch on its feasible vertices.
  ui target = 0;
  while (target < hitSets.size() && covered[target])
    target++;
  if (target == hitSets.size())
    return;

  for (ui v : hitSets[target]) {
    if (find(current.begin(), current.end(), v) != current.end())
      continue;
    bool inE = binary_search(E.begin(), E.end(), v);
    if (!inE)
      continue;

    bool connected = true;
    for (ui u : current) {
      if (!binary_search(adjList[u].begin(), adjList[u].end(), v)) {
        connected = false;
        break;
      }
    }
    if (!connected)
      continue;

    vector<ui> newlyCovered;
    for (ui i = 0; i < hitSets.size(); i++) {
      if (!covered[i] && binary_search(hitSets[i].begin(), hitSets[i].end(), v)) {
        covered[i] = 1;
        newlyCovered.push_back(i);
      }
    }

    current.push_back(v);
    backtrackSearch(E, hitSets, current, covered,
                    coveredCount + (ui)newlyCovered.size(), start + 1,
                    bestSize, solutions);
    current.pop_back();

    for (ui i : newlyCovered)
      covered[i] = 0;
  }
}

vector<vector<ui>>
ReorderSib::greedyApproximation(const vector<ui> &E,
                                const vector<vector<ui>> &hitSets) {
  vector<ui> S;
  vector<char> covered(hitSets.size(), 0);
  ui coveredCount = 0;

  while (coveredCount < hitSets.size()) {
    ui bestVertex = numeric_limits<ui>::max();
    ui bestGain = 0;

    for (ui v : E) {
      if (find(S.begin(), S.end(), v) != S.end())
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
        if (!covered[i] && binary_search(hitSets[i].begin(), hitSets[i].end(), v))
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

  ull fullMask = 0;
  for (ui i = 0; i < hitSets.size(); i++)
    fullMask |= (1ULL << i);

  vector<ull> masks(E.size(), 0);
  for (ui i = 0; i < E.size(); i++) {
    for (ui j = 0; j < hitSets.size(); j++) {
      if (binary_search(hitSets[j].begin(), hitSets[j].end(), E[i]))
        masks[i] |= (1ULL << j);
    }
  }

  vector<vector<ui>> solutions;
  vector<ui> current;

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
      for (ui v : current) {
        if (!binary_search(adjList[v].begin(), adjList[v].end(), E[i])) {
          connected = false;
          break;
        }
      }
      if (!connected)
        continue;

      current.push_back(E[i]);
      dfs(i + 1, remaining - 1, mask | masks[i]);
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

bool ReorderSib::branchSpaceInsideClique(const vector<ui> &M,
                                         const vector<ui> &E,
                                         const vector<ui> &C) {
  for (ui v : M) {
    if (!binary_search(C.begin(), C.end(), v))
      return false;
  }
  for (ui v : E) {
    if (!binary_search(C.begin(), C.end(), v))
      return false;
  }
  return true;
}

void ReorderSib::rCall(vector<vector<ui>> mustin,
                       vector<vector<ui>> expandTo, ui level,
                       vector<char> fullSkipCheck) {
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
    vector<ui> coveringCliques = collectCoveringCliques(baseM, level);
    bool hasCoveringCliques = !coveringCliques.empty();
    bool needsFullSkip = fullSkipCheck.empty() ? false : fullSkipCheck[0];
    vector<vector<ui>> siblingSets = generateSiblingSets(baseM, baseE, level);

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
          bool connectsToBase = true;
          for (ui u : baseMustin) {
            if (!binary_search(adjList[u].begin(), adjList[u].end(), v)) {
              connectsToBase = false;
              break;
            }
          }
          if (!connectsToBase)
            continue;

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

  for (ui i = 0; i < (ui)mustin.size(); i++) {
    vector<ui> R = mustin[i];
    vector<ui> Q = expandTo[i];
    bool done = false;
    enumerate(R, Q, mustin, expandTo, fullSkipCheck, i, level, done);
    if (done)
      break;
  }
}

void ReorderSib::enumerate(vector<ui> &R, vector<ui> &Q,
                           vector<vector<ui>> &mustin,
                           vector<vector<ui>> &expandTo,
                           vector<char> &fullSkipCheck, ui treeIndex,
                           ui level, bool &done) {
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

  if (Q.empty()) {
    if ((ui)R.size() > 2) {
      vector<ui> C = R;
      sort(C.begin(), C.end());
      bool alreadyFound =
          find(allCliques.begin(), allCliques.end(), C) != allCliques.end();
      if (!alreadyFound) {
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
      }

      done = true;
      vector<vector<ui>> newMustin;
      vector<vector<ui>> newExpandTo;
      vector<char> newFullSkipCheck;

      if (level == 0)
        cout << "here" << endl;
      else
        cout << "new here" << endl;

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
  cliqueCount = 0;
  maxCliqueSize = 0;
  checksCount = 0;
  allCliques.clear();
  foundLevel.clear();
  for (ui v = 0; v < n; v++)
    cliquesByVertex[v].clear();

  vector<vector<ui>> mustin;
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

  cout << "ReorderSib: cliques=" << cliqueCount
       << "  maxSize=" << maxCliqueSize << "  checks=" << checksCount
       << "  time=" << ms << " ms" << endl;
}
