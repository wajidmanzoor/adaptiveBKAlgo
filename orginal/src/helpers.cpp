#include "../inc/helpers.h"
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
  maxCliqueSize = 0; // Will be updated during search

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

bool PivotBK::isRedundantX(ui x, const vector<ui> &P) {
  // Check if x is connected to all vertices in P
  for (ui p : P) {
    if (!isConnected(x, p)) {
      return false;
    }
  }
  return true;
}
ui PivotBK::choosePivot(const vector<ui> &P, const vector<ui> &X) {
  ui bestPivot = P.empty() ? (X.empty() ? 0 : X[0]) : P[0];
  ui maxElimination = 0;

  // Check vertices in P to find the one that maximizes |P ∩ N(u)|
  for (ui u : P) {
    ui elimination = intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }

  // Check vertices in X to find the one that maximizes |P ∩ N(u)|
  for (ui u : X) {
    ui elimination = intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }

  // maximize |P ∩ N(u)|
  return bestPivot;
}
void PivotBK::applyDegreePruning(vector<ui> &P, const vector<ui> &R) {
  // Remove vertices with degree less than current clique size
  P.erase(remove_if(P.begin(), P.end(),
                    [this, &R](ui v) { return adjList[v].size() < R.size(); }),
          P.end());
}

void PivotBK::applyNeighborhoodPruning(vector<ui> &P, const vector<ui> &R) {
  // Remove vertices that can't form cliques larger than current max
  if (maxCliqueSize == 0)
    return;
  P.erase(remove_if(P.begin(), P.end(),
                    [this, &R, &P](ui v) {
                      ui potentialSize =
                          intersect(P, adjList[v]).size() + R.size();
                      return potentialSize < maxCliqueSize;
                    }),
          P.end());
}
void PivotBK::bronKerboschRecursive(vector<ui> &R, vector<ui> &P,
                                    vector<ui> &X) {
  // Basic pruning: check if P and X are empty
  if (isEmpty(P) && isEmpty(X)) {
    // Found a maximal clique
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());

    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui v : R) {
        cout << v << " ";
      }
      cout << "}" << endl;
    }
    return;
  }

  // Apply pruning rules
  // any vertx in P with degree less than current clique size (size of R) is
  // removed
  applyDegreePruning(P, R);

  // If P is empty after pruning, return
  if (isEmpty(P)) {
    return;
  }
  // Choose pivot using Tomita et al. strategy
  ui pivot = choosePivot(P, X);

  vector<ui> P_copy = P; // Copy to iterate over

  for (ui v : P_copy) {
    // Skip if v is a neighbor of the pivot
    // iter only to P \ N(pivot)
    if (isConnected(v, pivot)) {
      continue;
    }
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
void PivotBK::findAllMaximalCliques() {
  // Initialize sets
  vector<ui> R; // Current clique (empty)
  vector<ui> P; // All vertices as candidates
  vector<ui> X; // Excluded set (empty)

  // Fill P with all vertices
  for (ui i = 0; i < n; i++) {
    P.push_back(i);
  }

  cliqueCount = 0;
  maxCliqueSize = 0;
  bronKerboschRecursive(R, P, X);

  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}

// Reordering Bron-Kerbosch Implementation
ReorderBK::ReorderBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  visited.resize(n, false);

  // Fill adjacency lists from graph
  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < n) {
        adjList[i].push_back(neighbor);
      }
    }
    // Sort for efficient connectivity checks
    sort(adjList[i].begin(), adjList[i].end());
  }
}
bool ReorderBK::isConnected(ui u, ui v) const {
  // Binary search v in sorted adjacency list of u. I.e if edg (u, v) exists
  return binary_search(adjList[u].begin(), adjList[u].end(), v);
}
bool ReorderBK::canExtend(const vector<ui> &R, ui vertex) const {
  // Check if adding vertex to R maintains clique property
  for (ui v : R) {
    if (!isConnected(vertex, v)) {
      return false;
    }
  }
  return true;
}

void ReorderBK::findAllMaximalCliques() {
  cliqueCount = 0;
  maxCliqueSize = 0;

  cout << "Starting Reordering Bron-Kerbosch Algorithm..." << endl;

  vector<ui> ExpandFrom, ExpandMid, ExpandTo;

  ExpandFrom.push_back(0); // Start from vertex 0
  for (ui v = 1; v < n; v++) {
    ExpandTo.push_back(v);
    ExpandFrom.push_back(v);
  }

  rCall(ExpandFrom, ExpandMid, ExpandTo);

  cout << "\n=== Algorithm Complete ===" << endl;
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}

void ReorderBK::rCall(vector<ui> &ExpandFrom, vector<ui> &ExpandMid,
                      vector<ui> &ExpandTo) {
  if (debug) {
    cout << "\nRCall State:" << endl;
    cout << "  ExpandFrom: { ";
    for (ui v : ExpandFrom)
      cout << v << " ";
    cout << "}" << endl;

    cout << "  ExpandMid: { ";
    for (ui v : ExpandMid)
      cout << v << " ";
    cout << "}" << endl;

    cout << "  ExpandTo: { ";
    for (ui v : ExpandTo)
      cout << v << " ";
    cout << "}" << endl;
  }
  if (ExpandFrom.empty() || ExpandTo.empty()) {

    cout << " Empty From and To" << endl;
    return;
  }
  ui vertex = ExpandFrom[0];
  vector<ui> clique;
  clique.push_back(vertex);

  if (!ExpandMid.empty()) {
    cout << "Mid not empty" << endl;
    if (canExtend(clique, ExpandMid[0])) {
      clique.push_back(ExpandMid[0]);
      cout << "added " << ExpandMid[0] << " to p clique" << endl;
    } else {

      cout << "Expand using new Mid" << endl;
      vector<ui> newExpandMid;
      newExpandMid.push_back(ExpandTo[0]);
      ExpandTo.erase(ExpandTo.begin());
      rCall(ExpandFrom, newExpandMid, ExpandTo);
      return;
    }
  }

  for (ui v : ExpandTo) {
    if (canExtend(clique, v)) {
      clique.push_back(v);
      cout << "Added to clique " << v << endl;
    }
  }
  if (clique.size() > 2) {
    // Does that Mean what is left is a clique?
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)clique.size());
    if (debug) {
      cout << "Maximal Clique Found: { ";
      for (ui x : clique)
        cout << x << " ";
      cout << "}" << endl;
    }
    vector<ui> newExpandFrom, newExpandMid, newExpandTo;
    ExpandFrom.erase(ExpandFrom.begin());

    for (ui v : ExpandFrom) {
      if (!binary_search(clique.begin(), clique.end(), v) && !visited[v]) {
        newExpandFrom.push_back(v);
      }
    }
    ui ind = 0;
    while (ind < clique.size()) {
      if (visited[clique[ind]]) {
        ind++;
        continue;
      }
      newExpandMid.push_back(clique[ind]);
      ind++;
      break;
    }
    for (ui v = ind; v < clique.size(); v++) {
      if (!visited[clique[ind]])
        newExpandTo.push_back(clique[v]);
    }

    cout << "Expand after reorder" << endl;
    visited[vertex] = true;
    rCall(newExpandFrom, newExpandMid, newExpandTo);
    return;
  } else {
    if (!ExpandMid.empty()) {
      cout << "no clique found and Mid not empty expand" << endl;
      vector<ui> newExpandMid;
      newExpandMid.push_back(ExpandTo[0]);
      ExpandTo.erase(ExpandTo.begin());
      rCall(ExpandFrom, newExpandMid, ExpandTo);
      return;
    } else {
      cout << "No clique found and mid empty" << endl;
      vector<ui> newExpandFrom;
      newExpandFrom.push_back(ExpandTo[0]);
      ExpandTo.erase(ExpandTo.begin());
      rCall(newExpandFrom, ExpandMid, ExpandTo);
      return;
    }
  }
}