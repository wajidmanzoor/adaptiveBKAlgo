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

    if (R.size() > 2) {
      ofstream outfile("bk_pivot_maximal_cliques.txt", std::ios::app);
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

// Latest Reordering Bron-Kerbosch Implementation
ReorderBK2::ReorderBK2(Graph &g) {
  graph = g;
  n = g.n;
  adjList.resize(n);
  adjList2.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;

  // Fill adjacency lists from graph
  for (ui i = 0; i < n; i++) {
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if ((neighbor < n)) {
        // All neighbors
        adjList[i].push_back(neighbor);
      }
      if ((neighbor < n) && (neighbor > i)) {
        // Only greater than vertex neighbors
        adjList2[i].push_back(neighbor);
      }
    }
    // Sort for efficient connectivity checks
    sort(adjList[i].begin(), adjList[i].end());
    sort(adjList2[i].begin(), adjList2[i].end());
  }
}

vector<ui> ReorderBK2::intersect(vector<ui> vector1, vector<ui> vector2) {

  // Intersection between vector1 and vector2
  vector<char> seen(n, 0);

  for (ui x : vector1)
    seen[x] = 1;

  std::vector<ui> C;
  for (ui y : vector2)
    if (seen[y])
      C.push_back(y);
  return C;
}
vector<ui> ReorderBK2::unionSet(vector<ui> vector1, vector<ui> vector2) {

  // Union of vector1 and vector2
  vector<char> seen(n, 0);
  vector<ui> U;

  // Add elements from vector1
  for (ui x : vector1) {
    if (!seen[x]) {
      seen[x] = 1;
      U.push_back(x);
    }
  }

  // Add elements from vector2
  for (ui y : vector2) {
    if (!seen[y]) {
      seen[y] = 1;
      U.push_back(y);
    }
  }

  return U;
}

vector<ui> ReorderBK2::setDifference(vector<ui> A, vector<ui> B) {

  // Compute A - B
  vector<char> seen(n, 0);

  // Mark elements of B
  for (ui x : B)
    seen[x] = 1;

  vector<ui> C;
  for (ui y : A)
    if (!seen[y])
      C.push_back(y);

  return C;
}

void ReorderBK2::findAllMaximalCliques() {
  cliqueCount = 0;
  maxCliqueSize = 0;

  cout << "Starting Reordering Bron-Kerbosch Algorithm..." << endl;

  // Initialize heads and expandTo.
  vector<vector<ui>> mustin;
  vector<vector<ui>> expandTo;

  // Fill mustin with all vertices as individual trees and expandTo with their
  // greater than neighbors
  for (ui v = 0; v < n; v++) {
    vector<ui> temp;
    temp.push_back(v);
    mustin.push_back(temp);
    // TODO: no need as already in adjList2 with only greater vertices in
    // constructor replace with adjList2
    auto it = std::upper_bound(adjList[v].begin(), adjList[v].end(), v);
    expandTo.emplace_back(it, adjList[v].end());
  }

  // Start the recursive search.
  rCall(mustin, expandTo, 0, 0);

  cout << "\n=== Algorithm Complete ===" << endl;
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}

void ReorderBK2::rCall(vector<vector<ui>> &mustin, vector<vector<ui>> &expandTo,
                       ui level, ui enlevel) {

  cout << "RCall Level: " << level << endl;
  cout << "\nRCall State:" << endl;
  cout << "   Mustin Trees: " << endl;
  for (vector<ui> tree : mustin) {
    cout << "     { ";
    for (ui v : tree) {
      cout << v << " ";
    }
    cout << "}" << endl;
  }
  cout << endl << "   ExpandTo: " << endl;
  for (vector<ui> tree : expandTo) {
    cout << "     { ";
    for (ui v : tree) {
      cout << v << " ";
    }
    cout << "}" << endl;
  }
  cout << endl;

  if (mustin.empty()) {
    return;
  }

  // Partial clique candidates
  vector<ui> Q;

  // Partial clique
  vector<ui> R;

  // Falg that represents if we find a maximal clique in this rCall?
  // If yes, we do not need to continue further in this rCall, and trees will be
  // updated and new call will be made. If no, we continue to enemurate trees in
  // this rCall
  bool moveToNext = true;

  // Index to track which tree we are processing
  ui index = 0;
  for (ui i = 0; i < mustin.size(); i++) {

    // Stop if maximal clique already found
    if (!moveToNext)
      return;

    // Skip empty trees
    if (mustin[i].empty())
      continue;
    // Initialize R and Q for this tree
    R.clear();
    for (ui v : mustin[i]) {
      R.push_back(v);
    }
    Q.clear();
    Q = expandTo[i];

    // Flag to indicate if maximal clique found in enemurate. Usefull to stop
    // further enemuration in this rCall
    bool flag = false;
    enemurate(R, Q, mustin, expandTo, moveToNext, flag, index, level, enlevel);
    index++;
  }
  return;
}

void ReorderBK2::enemurate(vector<ui> &R, vector<ui> &Q,
                           vector<vector<ui>> &mustin,
                           vector<vector<ui>> &expandTo, bool &moveToNext,
                           bool &flag, ui index, ui level, ui enlevel) {

  cout << " rcall " << level << "  Enumerate Level: " << enlevel << endl;
  cout << "     Enemurate Call:" << endl;
  cout << "       Current R: ";
  for (ui v : R) {
    cout << v << " ";
  }
  cout << endl;
  cout << "       Current Q: ";
  for (ui v : Q) {
    cout << v << " ";
  }
  cout << endl;
  cout << "       MoveToNext: " << moveToNext << endl;
  cout << "       Current Index: " << index << endl;

  // No more enumeration left, Maximal clique found
  if (Q.empty()) {
    // Maximal clique found is not a edge or vertex
    if (R.size() > 2) {

      // set flag to true to stop further enemuration in this rCall
      flag = true;

      // Update clique counts
      cliqueCount++;
      maxCliqueSize = max(maxCliqueSize, R.size());

      // Set moveToNext to false to stop further rCall (tree expansions) as tree
      // need to be reordered
      moveToNext = false;
      if (debug) {
        cout << endl << "************************" << endl;
        cout << "Maximal Clique: { ";
        for (ui v : R) {
          cout << v << " ";
        }
        cout << "}";
        cout << endl << "************************" << endl;
      }

      // Reordering logic to be added here

      // Create new trees (mustin and expandTo) for next rCall
      vector<vector<ui>> newMustin;
      vector<vector<ui>> newExpandTo;

      // removes trees who mustin vertices are all in maximal Clique(R).
      // For remaining trees, update their expandTo sets by including Verticies
      // in of maximal Clique (R), if they are connected to all vertices in
      // mustin.
      for (ui i = index; i < mustin.size(); i++) {
        vector<ui> tree = mustin[i];
        ui len = tree.size();

        // Check how many vertices of tree are in R
        for (ui v : tree) {
          if (find(R.begin(), R.end(), v) < R.end()) {
            len--;
          }
        }
        // If some vertices remain, keep this tree for further expansion
        if (len > 0) {
          newMustin.push_back(tree);
          vector<ui> temp2 = R;
          temp2 = unionSet(temp2, expandTo[i]);
          for (ui v : tree) {
            temp2 = intersect(temp2, adjList[v]);
          }
          newExpandTo.push_back(temp2);
        }
      }

      cout << "       New Mustin Trees after removing clique vertices: "
           << endl;
      for (vector<ui> tree : newMustin) {
        cout << "         { ";
        for (ui v : tree) {
          cout << v << " ";
        }
        cout << "}" << endl;
      }
      cout << endl;

      cout << "       New ExpandTo after removing clique vertices: " << endl;
      for (vector<ui> tree : newExpandTo) {
        cout << "         { ";
        for (ui v : tree) {
          cout << v << " ";
        }
        cout << "}" << endl;
      }
      cout << endl;

      // For each tree, expand mustin by adding each vertex from expandTo one
      // at a time and updating expandTo accordingly. This dicides each tree
      // into its individual subtrees for next rCall.
      for (ui i = 0; i < newMustin.size(); i++) {
        vector<vector<ui>> mustin_;
        vector<vector<ui>> expand_;
        if (newExpandTo[i].empty()) {
          mustin_.push_back(newMustin[i]);
          expand_.push_back(newExpandTo[i]);
        }

        for (ui v : newExpandTo[i]) {
          vector<ui> temp = newMustin[i];
          temp.push_back(v);
          mustin_.push_back(temp);
          expand_.push_back(intersect(newExpandTo[i], adjList2[v]));
        }

        rCall(mustin_, expand_, level + 1, 0);
      }

      return;
    } else {
      return; // clique too small, edge or vertex: move on
    }
  }

  cout << "       Flag at start of enemurate: " << flag << endl;
  for (ui v : Q) {

    // Stop if maximal clique already found
    if (flag)
      return;
    // Add v to R: Partial Clique
    R.push_back(v);

    vector<ui> Q_;
    // Q_ = Q ∩ N(v): New possible candidates i.e. verticies connected to all
    // verticies in Partial clique R. Use adjList2 for greater than vertices
    // neighbor only
    Q_ = intersect(Q, adjList2[v]);

    // Recursive call to enemurate with expanded R and new Q_
    enemurate(R, Q_, mustin, expandTo, moveToNext, flag, index, level,
              enlevel + 1);

    // Backtrack: remove v from R
    R.pop_back();
  }
}

/// new implementation of ReorderBK2 class
// Reordering Bron-Kerbosch Implementation

ReorderBK::ReorderBK(Graph &g) {
  graph = g;
  n = g.n;
  adjList.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;
  visited.resize(n, false);
  status.resize(n, false);

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

vector<ui> ReorderBK::intersect(vector<ui> vector1, vector<ui> vector2) {

  // Intersection between vector1 and vector2
  vector<char> seen(n, 0);

  for (ui x : vector1)
    seen[x] = 1;

  std::vector<ui> C;
  for (ui y : vector2)
    if (seen[y])
      C.push_back(y);
  return C;
}

void ReorderBK::enemurate(vector<ui> &R, vector<ui> &Q, vector<ui> &heads,
                          vector<ui> &expandTo) {

  // If Already visited, return
  if (visited[R[0]]) {
    // redendant call
    return;
  }
  cout << "   Checking node: " << endl;
  cout << "   Current R: ";

  for (ui x : R) {
    cout << x << " ";
  }
  cout << endl << "   Current Q: ";
  for (ui x : Q) {
    cout << x << " ";
  }
  cout << endl;

  // If Q is empty i.e if no more enumeration left, check R for clique and
  // reorder or return
  if (Q.empty()) {

    // If maximal clique size > 2, process clique
    if (R.size() > 2) {

      cliqueCount++;
      maxCliqueSize = max(maxCliqueSize, (ui)R.size());

      if (debug) {
        cout << endl << "************************" << endl;
        cout << "Maximal Clique: { ";
        for (ui v : R) {
          cout << v << " ";
        }
        cout << "}";
        cout << endl << "************************" << endl;
      }

      // If not the last Head, perform reorder call
      if (heads.size() > 1) {

        // Reordering

        // Prepare new heads
        vector<ui> newHeads;
        heads.erase(heads.begin());

        // Update status and visited for vertices in found clique
        for (ui v : R) {
          status[v] = cliqueCount;
        }
        for (ui v : R) {
          visited[v] = true;
        }

        // Fill new heads excluding clique vertices
        for (ui v : heads) {
          if (status[v] != cliqueCount) {
            newHeads.push_back(v);
          }
        }

        // Expand with new head
        rCall(newHeads, expandTo);

        // clique found and reorder call
        return;
      } else {
        // No reorder if we are expanding using the last head.
        if (debug) {
          cout << endl << "************************" << endl;
          cout << "Maximal Clique: { ";
          for (ui v : R) {
            cout << v << " ";
          }
          cout << "}";
          cout << endl << "************************" << endl;
        }
        return; // head fully exhausted as no need for redorder as last head
                // is explored
      }

    } else {

      return; // clique too small, edge or vertex: move on
    }
  }

  // Continue enemuration by expanding each vertex in Q
  for (ui v : Q) {

    // Add v to R: Partial Clique
    R.push_back(v);

    vector<ui> Q_;
    // Q_ = Q ∩ N(v): New possible candidates i.e. verticies connected to all
    // verticies in Partial clique R
    Q_ = intersect(Q, adjList[v]);

    // Recursive call to enemurate with expanded R and new Q_
    enemurate(R, Q_, heads, expandTo);

    // Backtrack: remove v from R
    R.pop_back();
  }
}

void ReorderBK::rCall(vector<ui> heads, vector<ui> expandTo) {
  if (debug) {
    cout << "\nRCall State:" << endl;
    cout << "head : { ";
    for (ui v : heads)
      cout << v << " ";
    cout << "}" << endl;

    cout << "Status: { ";
    for (ui v : status)
      cout << v << " ";
    cout << "}" << endl;

    cout << "ExpandTo: { ";
    for (ui v : expandTo)
      cout << v << " ";
    cout << "}" << endl;

    cout << "Visited: ";
    for (bool val : visited) {
      cout << val << " ";
    }
    cout << endl;
  }

  // Base case: no more heads to process
  if (heads.empty()) {
    return;
  }

  // Process the one head at a time
  ui vertex = heads[0];

  // Initialize R (partial clique) and Q (Possible candidates) for enemuration
  vector<ui> R;
  R.push_back(vertex);
  vector<ui> Q;

  // Q = expandTo ∩ N(vertex)
  Q = intersect(adjList[vertex], expandTo);

  // Start enemuration with current head: Either finds clique and reorder or
  // marks visited and move to next head
  enemurate(R, Q, heads, expandTo);

  // If no clique found by expanding this head, mark as visited and proceed
  if (!visited[vertex]) {
    cout << "vertex " << vertex << " visited no clique " << endl;

    // Move to next head
    heads.erase(heads.begin());
    // Mark vertex as visited
    visited[vertex] = true;

    // remove vertex from expandTo
    vector<ui> newExpandTo;
    for (ui v : expandTo) {
      if (v != vertex) {
        newExpandTo.push_back(v);
      }
    }

    // Expand new heads with updated expandTo
    rCall(heads, newExpandTo);
  }
}

void ReorderBK::findAllMaximalCliques() {
  cliqueCount = 0;
  maxCliqueSize = 0;

  cout << "Starting Reordering Bron-Kerbosch Algorithm..." << endl;

  // Initialize heads and expandTo.
  vector<ui> heads;
  vector<ui> expandTo;

  // Fill heads and expandTo with all vertices
  for (ui v = 0; v < n; v++) {
    expandTo.push_back(v);
    heads.push_back(v);
  }

  // Start the recursive search one head at a time
  rCall(heads, expandTo);

  cout << "\n=== Algorithm Complete ===" << endl;
  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;
}
