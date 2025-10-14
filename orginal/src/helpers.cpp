#include "../inc/helpers.h"

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

vector<ui> PivotBK::intersect(const vector<ui> &set1, const vector<ui> &neighbors) {
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

bool PivotBK::isEmpty(const vector<ui> &set) { 
  return set.empty(); 
}

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
  
  // Check vertices in P
  for (ui u : P) {
    ui elimination = intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }
  
  // Check vertices in X
  for (ui u : X) {
    ui elimination = intersect(P, adjList[u]).size();
    if (elimination > maxElimination) {
      maxElimination = elimination;
      bestPivot = u;
    }
  }
  
  return bestPivot;
}

void PivotBK::applyDegreePruning(vector<ui> &P, const vector<ui> &R) {
  // Remove vertices with degree less than current clique size
  P.erase(remove_if(P.begin(), P.end(), 
    [this, &R](ui v) { 
      return adjList[v].size() < R.size(); 
    }), P.end());
}

void PivotBK::applyNeighborhoodPruning(vector<ui> &P, const vector<ui> &R) {
  // Only apply if we have found at least one clique
  if (maxCliqueSize == 0) return;
  
  // Remove vertices that can't form cliques larger than current max
  P.erase(remove_if(P.begin(), P.end(), 
    [this, &R, &P](ui v) { 
      ui potentialSize = intersect(P, adjList[v]).size() + R.size();
      return potentialSize < maxCliqueSize; 
    }), P.end());
}

void PivotBK::bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X) {
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

  // Basic pruning: remove redundant vertices from X
  // X.erase(remove_if(X.begin(), X.end(), 
  //   [this, &P](ui x) { return isRedundantX(x, P); }), X.end());

  // Apply pruning rules
  applyDegreePruning(P, R);
  // applyNeighborhoodPruning(P, R); // Disabled - too aggressive

  // If P is empty after pruning, return
  if (isEmpty(P)) {
    return;
  }

  // Choose pivot using Tomita et al. strategy
  ui pivot = choosePivot(P, X);
  
  // Process each vertex in P that is NOT a neighbor of the pivot
  // This is the key optimization: we don't need to process neighbors of pivot
  // because they will be processed in other branches
  vector<ui> P_copy = P; // Copy to iterate over
  for (ui v : P_copy) {
    // Skip if v is a neighbor of the pivot
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

// ReorderBK Implementation
ReorderBK::ReorderBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  foundCliques.clear();

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

vector<ui> ReorderBK::intersect(const vector<ui> &set1, const vector<ui> &neighbors) {
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

bool ReorderBK::isEmpty(const vector<ui> &set) { 
  return set.empty(); 
}

set<ui> ReorderBK::vectorToSet(const vector<ui> &vec) {
  return set<ui>(vec.begin(), vec.end());
}

bool ReorderBK::isSubset(const set<ui> &subset, const set<ui> &superset) {
  // Check if subset is a subset of superset
  return includes(superset.begin(), superset.end(), subset.begin(), subset.end());
}

bool ReorderBK::isSubsetOfFoundClique(const vector<ui> &candidate) {
  set<ui> candidateSet = vectorToSet(candidate);
  
  for (const auto& foundClique : foundCliques) {
    if (isSubset(candidateSet, foundClique)) {
      if (debug) {
        cout << "Found subset: { ";
        for (ui v : candidate) cout << v << " ";
        cout << "} is subset of { ";
        for (ui v : foundClique) cout << v << " ";
        cout << "}" << endl;
      }
      return true;
    }
  }
  return false;
}

void ReorderBK::bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X) {
  // Note: We removed early subset checking to avoid being too aggressive
  // Subset checking is now only done when we find a maximal clique

  // Basic pruning: check if P and X are empty
  if (isEmpty(P) && isEmpty(X)) {
    // Check if this clique is a subset of any found clique
    if (!foundCliques.empty() && isSubsetOfFoundClique(R)) {
      if (debug) {
        cout << "Skipping duplicate clique: { ";
        for (ui v : R) cout << v << " ";
        cout << "}" << endl;
      }
      return; // Skip this clique - it's a subset of a found one
    }
    
    // Found a maximal clique
    cliqueCount++;
    
    // Add to found cliques set
    set<ui> newClique = vectorToSet(R);
    foundCliques.insert(newClique);
    
    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui v : R) {
        cout << v << " ";
      }
      cout << "}" << endl;
    }
    return;
  }

  // Process each vertex in P
  vector<ui> P_copy = P; // Copy to iterate over
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

void ReorderBK::findAllMaximalCliques() {
  // Initialize sets
  vector<ui> R; // Current clique (empty)
  vector<ui> P; // All vertices as candidates
  vector<ui> X; // Excluded set (empty)

  // Fill P with all vertices
  for (ui i = 0; i < n; i++) {
    P.push_back(i);
  }

  cliqueCount = 0;
  foundCliques.clear();
  bronKerboschRecursive(R, P, X);

  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
}

// PureReorderBK Implementation
PureReorderBK::PureReorderBK(Graph &g) : G(g.n) {
  // Convert adjacency list to bitset representation
  for (ui i = 0; i < g.n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < g.n) {
        G.addEdge(i, neighbor);
      }
    }
  }
  
  // Initialize global order (all vertices)
  globalOrder.clear();
  for (ui i = 0; i < g.n; i++) {
    globalOrder.push_back(i);
  }
  
  cliqueCount = 0;
  foundCliques.clear();
}

bool PureReorderBK::canAdd(const vector<ui> &R, ui v) {
  // Check if v is connected to all vertices in R
  for (ui u : R) {
    if (!G.connected(u, v)) {
      return false;
    }
  }
  return true;
}

void PureReorderBK::reorderAfterClique(const vector<ui> &R) {
  // Create a vector to mark which vertices are in the found clique
  vector<bool> inR(G.n, false);
  for (ui v : R) {
    inR[v] = true;
  }
  
  // Separate vertices into front (not in clique) and back (in clique)
  vector<ui> front, back;
  front.reserve(G.n);
  back.reserve(G.n);
  
  for (ui v : globalOrder) {
    if (inR[v]) {
      back.push_back(v);
    } else {
      front.push_back(v);
    }
  }
  
  // Rebuild global order: front vertices first, then clique vertices
  globalOrder.clear();
  globalOrder.insert(globalOrder.end(), front.begin(), front.end());
  globalOrder.insert(globalOrder.end(), back.begin(), back.end());
  
  if (debug) {
    cout << "Reordered: [";
    for (ui v : globalOrder) cout << v << " ";
    cout << "]" << endl;
  }
}

set<ui> PureReorderBK::vectorToSet(const vector<ui> &vec) {
  return set<ui>(vec.begin(), vec.end());
}

bool PureReorderBK::isSubsetOfFoundClique(const vector<ui> &candidate) {
  set<ui> candidateSet = vectorToSet(candidate);
  
  for (const auto& foundClique : foundCliques) {
    if (includes(foundClique.begin(), foundClique.end(), candidateSet.begin(), candidateSet.end())) {
      return true;
    }
  }
  return false;
}

void PureReorderBK::reportClique(const vector<ui> &R) {
  // Check if this clique is a subset of any found clique
  if (!foundCliques.empty() && isSubsetOfFoundClique(R)) {
    if (debug) {
      cout << "Skipping duplicate clique: { ";
      for (ui v : R) cout << v << " ";
      cout << "}" << endl;
    }
    return; // Skip this clique - it's a subset of a found one
  }
  
  cliqueCount++;
  
  // Add to found cliques set
  set<ui> newClique = vectorToSet(R);
  foundCliques.insert(newClique);
  
  if (debug) {
    cout << "Maximal Clique: { ";
    for (ui v : R) {
      cout << v << " ";
    }
    cout << "}" << endl;
  }
}

void PureReorderBK::enumerate(vector<ui> &R, ui startIdx) {
  bool expanded = false;
  
  // Try to add any vertex from current order starting at startIdx
  for (ui i = startIdx; i < globalOrder.size(); i++) {
    ui v = globalOrder[i];
    if (!canAdd(R, v)) continue;
    
    // Extend clique with vertex v
    R.push_back(v);
    enumerate(R, i + 1);  // Continue from i+1
    R.pop_back();
    expanded = true;
  }
  
  // If we couldn't add anything, R is maximal (unless empty)
  if (!expanded && !R.empty()) {
    reportClique(R);
    // Reorder after finding a clique
    reorderAfterClique(R);
  }
}

void PureReorderBK::findAllMaximalCliques() {
  cliqueCount = 0;
  
  // Start with empty clique
  vector<ui> R;
  R.reserve(G.n);
  
  // Begin enumeration from index 0
  enumerate(R, 0);
  
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
}