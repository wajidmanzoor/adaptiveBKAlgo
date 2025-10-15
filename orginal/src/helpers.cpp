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

// AdaptiveSkipPivotBK Implementation
AdaptiveSkipPivotBK::AdaptiveSkipPivotBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  
  // Initialize skip mask
  skip_mask.resize(n, false);
  
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
  
  // Initialize global order - start with degree-based ordering
  global_order.resize(n);
  iota(global_order.begin(), global_order.end(), 0);
  
  // Sort by degree (descending) for better initial ordering
  sort(global_order.begin(), global_order.end(), 
       [this](ui a, ui b) { return adjList[a].size() > adjList[b].size(); });
}

vector<ui> AdaptiveSkipPivotBK::intersect(const vector<ui> &set1, const vector<ui> &neighbors) {
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

bool AdaptiveSkipPivotBK::isEmpty(const vector<ui> &set) { 
  return set.empty(); 
}

bool AdaptiveSkipPivotBK::isConnected(ui u, ui v) {
  return binary_search(adjList[u].begin(), adjList[u].end(), v);
}

ui AdaptiveSkipPivotBK::choosePivot(const vector<ui> &P, const vector<ui> &X) {
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

void AdaptiveSkipPivotBK::applyDegreePruning(vector<ui> &P, const vector<ui> &R) {
  // Remove vertices with degree less than current clique size
  P.erase(remove_if(P.begin(), P.end(),
                    [this, &R](ui v) { return adjList[v].size() < R.size(); }),
          P.end());
}

bool AdaptiveSkipPivotBK::isClique(const vector<ui> &R) const {
  for (size_t i = 0; i < R.size(); i++) {
    for (size_t j = i + 1; j < R.size(); j++) {
      if (!binary_search(adjList[R[i]].begin(), adjList[R[i]].end(), R[j])) {
        return false;
      }
    }
  }
  return true;
}

void AdaptiveSkipPivotBK::reorderAfterClique(const vector<ui> &clique) {
  // For now, let's focus on reordering without aggressive skipping
  // Mark clique vertices as processed (but don't skip them completely)
  for (ui v : clique) {
    skip_mask[v] = true;
  }
  
  // Reorder: unprocessed vertices first, processed vertices last
  vector<ui> unprocessed, processed;
  for (ui v : global_order) {
    if (skip_mask[v]) {
      processed.push_back(v);
    } else {
      unprocessed.push_back(v);
    }
  }
  
  // Update global order to prioritize unprocessed vertices
  global_order.clear();
  global_order.insert(global_order.end(), unprocessed.begin(), unprocessed.end());
  global_order.insert(global_order.end(), processed.begin(), processed.end());
  
  if (debug) {
    cout << "Reordered vertices: unprocessed=[";
    for (ui v : unprocessed) cout << v << " ";
    cout << "], processed=[";
    for (ui v : processed) cout << v << " ";
    cout << "]" << endl;
  }
}

vector<ui> AdaptiveSkipPivotBK::getOrderedCandidates(const vector<ui> &P, ui start_idx) const {
  vector<ui> ordered_candidates;
  
  // Follow global order starting from start_idx, skip covered vertices
  for (ui idx = start_idx; idx < global_order.size(); idx++) {
    ui v = global_order[idx];
    if (skip_mask[v]) continue;  // Skip covered vertices
    
    // Check if vertex is in P
    if (find(P.begin(), P.end(), v) != P.end()) {
      ordered_candidates.push_back(v);
    }
  }
  
  return ordered_candidates;
}

void AdaptiveSkipPivotBK::adaptiveBronKerbosch(vector<ui> &R, vector<ui> &P, vector<ui> &X, ui start_idx) {
  // Base case: maximal clique found
  if (isEmpty(P) && isEmpty(X)) {
    cliqueCount++;
    maxCliqueSize = max(maxCliqueSize, (ui)R.size());
    
    if (debug) {
      cout << "Maximal Clique: { ";
      for (ui v : R) cout << v << " ";
      cout << "}" << endl;
    }
    
    // Apply skip-mask reordering strategy
    reorderAfterClique(R);
    return;
  }
  
  // Get candidates in global order - don't skip covered vertices in P construction
  // The skip happens during iteration to allow for proper P/X management
  vector<ui> ordered_candidates;
  for (ui idx = start_idx; idx < global_order.size(); idx++) {
    ui v = global_order[idx];
    if (find(P.begin(), P.end(), v) != P.end()) {
      ordered_candidates.push_back(v);
    }
  }
  
  // Apply degree pruning
  applyDegreePruning(ordered_candidates, R);
  if (ordered_candidates.empty()) return;
  
  // Choose pivot from ordered candidates + X
  ui pivot = choosePivot(ordered_candidates, X);
  
  // Process candidates in global order
  vector<ui> candidates_copy = ordered_candidates;
  for (ui v : candidates_copy) {
    // Skip if v is a neighbor of the pivot (standard pivoting)
    if (isConnected(v, pivot)) continue;
    
    // Skip if vertex is already covered by a maximal clique (skip-mask optimization)
    if (skip_mask[v]) continue;
    
    vector<ui> new_R = R;
    new_R.push_back(v);
    
    // Verify that new_R forms a valid clique (debug)
    if (!isClique(new_R)) {
      if (debug) {
        cout << "Skipping invalid clique: { ";
        for (ui x : new_R) cout << x << " ";
        cout << "}" << endl;
      }
      continue;
    }
    
    // Create P ∩ N(v) and X ∩ N(v)
    vector<ui> new_P = intersect(P, adjList[v]);
    vector<ui> new_X = intersect(X, adjList[v]);
    
    // Find next start index in global order
    ui next_start = find(global_order.begin(), global_order.end(), v) - global_order.begin() + 1;
    
    // Recursive call
    adaptiveBronKerbosch(new_R, new_P, new_X, next_start);
    
    // Move v from P to X
    P.erase(find(P.begin(), P.end(), v));
    X.push_back(v);
  }
}

void AdaptiveSkipPivotBK::findAllMaximalCliques() {
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
  
  cout << "Starting Adaptive Skip-Mask Bron-Kerbosch..." << endl;
  adaptiveBronKerbosch(R, P, X, 0);  // Start from index 0
  
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}

// SimpleAdaptiveBK Implementation - Pure enumeration with skip-mask reordering
SimpleAdaptiveBK::SimpleAdaptiveBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  
  // Initialize skip mask
  skip_mask.resize(n, false);
  
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
  
  // Initialize global order - simple sequential order
  global_order.resize(n);
  iota(global_order.begin(), global_order.end(), 0);
}

bool SimpleAdaptiveBK::isConnected(ui u, ui v) const {
  return binary_search(adjList[u].begin(), adjList[u].end(), v);
}

bool SimpleAdaptiveBK::isClique(const vector<ui> &R) const {
  for (size_t i = 0; i < R.size(); i++) {
    for (size_t j = i + 1; j < R.size(); j++) {
      if (!isConnected(R[i], R[j])) {
        return false;
      }
    }
  }
  return true;
}

void SimpleAdaptiveBK::handleClique(const vector<ui> &R) {
  cliqueCount++;
  maxCliqueSize = max(maxCliqueSize, (ui)R.size());
  
  if (debug) {
    cout << "Found maximal clique: { ";
    for (ui v : R) cout << v << " ";
    cout << "}" << endl;
  }
  
  // Mark covered vertices - but be more conservative
  // Only mark vertices that are definitely covered by this maximal clique
  for (ui v : R) {
    skip_mask[v] = true;
  }
  
  // Reorder: uncovered first, covered last
  vector<ui> uncovered, covered;
  for (ui v : global_order) {
    if (skip_mask[v]) {
      covered.push_back(v);
    } else {
      uncovered.push_back(v);
    }
  }
  
  global_order.clear();
  global_order.insert(global_order.end(), uncovered.begin(), uncovered.end());
  global_order.insert(global_order.end(), covered.begin(), covered.end());
  
  if (debug) {
    cout << "Reordered: uncovered=[";
    for (ui v : uncovered) cout << v << " ";
    cout << "], covered=[";
    for (ui v : covered) cout << v << " ";
    cout << "]" << endl;
  }
}

void SimpleAdaptiveBK::enumerate(vector<ui> &R, ui start_idx) {
  bool expanded = false;
  
  for (ui idx = start_idx; idx < global_order.size(); idx++) {
    ui v = global_order[idx];
    
    // Skip redundant vertices (already covered by maximal cliques)
    if (skip_mask[v]) continue;
    
    vector<ui> R_next = R;
    R_next.push_back(v);
    
    // Check if R_next forms a clique
    if (isClique(R_next)) {
      expanded = true;
      enumerate(R_next, idx + 1);
    }
  }
  
  // Leaf node → maximal clique found (cannot be extended further)
  if (!expanded && !R.empty()) {
    handleClique(R);
  }
}

void SimpleAdaptiveBK::findAllMaximalCliques() {
  vector<ui> R; // Start with empty clique
  
  cliqueCount = 0;
  maxCliqueSize = 0;
  
  cout << "Starting Simple Adaptive Enumeration..." << endl;
  enumerate(R, 0);  // Start from index 0
  
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}