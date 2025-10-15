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

// DepthFirstReorderBK Implementation
DepthFirstReorderBK::DepthFirstReorderBK(Graph &g) {
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  
  // Initialize algorithm state
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
  
  // Initialize global order - simple sequential order initially
  global_order.resize(n);
  iota(global_order.begin(), global_order.end(), 0);
}

bool DepthFirstReorderBK::isConnected(ui u, ui v) const {
  return binary_search(adjList[u].begin(), adjList[u].end(), v);
}

bool DepthFirstReorderBK::isClique(const vector<ui> &R) const {
  for (size_t i = 0; i < R.size(); i++) {
    for (size_t j = i + 1; j < R.size(); j++) {
      if (!isConnected(R[i], R[j])) {
        return false;
      }
    }
  }
  return true;
}

bool DepthFirstReorderBK::canExtend(const vector<ui> &R, ui vertex) const {
  // Check if adding vertex to R maintains clique property
  for (ui v : R) {
    if (!isConnected(vertex, v)) {
      return false;
    }
  }
  return true;
}

vector<ui> DepthFirstReorderBK::depthFirstExpand(ui start_vertex) {
  vector<ui> current_clique = {start_vertex};
  
  if (debug) {
    cout << "Starting depth-first expansion from vertex " << start_vertex << endl;
  }
  
  bool expanded = true;
  while (expanded) {
    expanded = false;
    
    // Try to extend current clique with ANY vertex in global order
    // This is the original logic - try all possible extensions
    for (ui v : global_order) {
      // Skip if vertex is already in current clique
      if (find(current_clique.begin(), current_clique.end(), v) != current_clique.end()) {
        continue;
      }
      
      // Check if adding v maintains clique property
      if (canExtend(current_clique, v)) {
        current_clique.push_back(v);
        expanded = true;
        
        if (debug) {
          cout << "Extended clique to: { ";
          for (ui x : current_clique) cout << x << " ";
          cout << "}" << endl;
        }
        break; // Continue with depth-first expansion
      }
    }
  }
  
  if (debug) {
    cout << "Found maximal clique: { ";
    for (ui v : current_clique) cout << v << " ";
    cout << "}" << endl;
  }
  
  return current_clique;
}

ui DepthFirstReorderBK::getNextStartingVertex() {
  // Find first unvisited vertex in global order
  for (ui v : global_order) {
    if (!visited[v]) {
      return v;
    }
  }
  return n; // All vertices visited
}

void DepthFirstReorderBK::reorderVertices() {
  vector<ui> tier1, tier2, tier3, tier4;
  
  // Classify vertices into 4 tiers
  for (ui v = 0; v < n; v++) {
    bool is_visited = visited[v];
    bool in_clique = false;
    
    // Check if vertex is in any found clique
    for (const auto& clique : found_cliques) {
      if (clique.find(v) != clique.end()) {
        in_clique = true;
        break;
      }
    }
    
    // Assign to appropriate tier
    if (!is_visited && !in_clique) {
      tier1.push_back(v); // Highest priority
    } else if (!is_visited && in_clique) {
      tier2.push_back(v);
    } else if (is_visited && !in_clique) {
      tier3.push_back(v);
    } else { // is_visited && in_clique
      tier4.push_back(v); // Lowest priority
    }
  }
  
  // Rebuild global order with new priorities
  global_order.clear();
  global_order.insert(global_order.end(), tier1.begin(), tier1.end());
  global_order.insert(global_order.end(), tier2.begin(), tier2.end());
  global_order.insert(global_order.end(), tier3.begin(), tier3.end());
  global_order.insert(global_order.end(), tier4.begin(), tier4.end());
  
  if (debug) {
    cout << "Reordered vertices:" << endl;
    cout << "  Tier 1 (not_visited, not_in_clique): [";
    for (ui v : tier1) cout << v << " ";
    cout << "]" << endl;
    cout << "  Tier 2 (not_visited, in_clique): [";
    for (ui v : tier2) cout << v << " ";
    cout << "]" << endl;
    cout << "  Tier 3 (visited, not_in_clique): [";
    for (ui v : tier3) cout << v << " ";
    cout << "]" << endl;
    cout << "  Tier 4 (visited, in_clique): [";
    for (ui v : tier4) cout << v << " ";
    cout << "]" << endl;
  }
}

void DepthFirstReorderBK::findAllMaximalCliques() {
  cliqueCount = 0;
  maxCliqueSize = 0;
  
  cout << "Starting Depth-First Reordering Algorithm..." << endl;
  
  // Main algorithm loop
  while (true) {
    ui start_vertex = getNextStartingVertex();
    
    // Termination condition: all vertices visited
    if (start_vertex == n) {
      break;
    }
    
    if (debug) {
      cout << "\n=== Phase " << (cliqueCount + 1) << ": Starting with vertex " << start_vertex << " ===" << endl;
    }
    
    // Mark vertex as visited
    visited[start_vertex] = true;
    
    // Perform depth-first expansion
    vector<ui> maximal_clique = depthFirstExpand(start_vertex);
    
    // Store the found maximal clique (only if it's new)
    set<ui> clique_set(maximal_clique.begin(), maximal_clique.end());
    
    // Check if this clique is already found
    if (found_cliques.find(clique_set) == found_cliques.end()) {
      found_cliques.insert(clique_set);
      cliqueCount++;
      maxCliqueSize = max(maxCliqueSize, (ui)maximal_clique.size());
      
      cout << "NEW maximal clique found: { ";
      for (ui v : maximal_clique) cout << v << " ";
      cout << "}" << endl;
    } else {
      if (debug) {
        cout << "Duplicate clique skipped: { ";
        for (ui v : maximal_clique) cout << v << " ";
        cout << "}" << endl;
      }
    }
    
    // Reorder vertices based on new information
    reorderVertices();
  }
  
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}