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
vector<ui> ReorderBK2::applyPreviousCliqueEffects(const vector<ui> &tree,
                                                  const vector<ui> &Q) {
  vector<ui> filtered = Q;

  // For every previously found maximal clique
  for (const auto &C : foundCliques) {
    bool affectsTree = false;

    // If this clique shares at least one vertex with current tree,
    // then it can constrain this tree
    for (ui x : tree) {
      if (find(C.begin(), C.end(), x) != C.end()) {
        affectsTree = true;
        break;
      }
    }

    if (!affectsTree)
      continue;

    // Remove clique vertices that are not already in the tree.
    // Those cannot be reused in expandTo for this tree.
    vector<ui> next;
    for (ui v : filtered) {
      bool inClique = (find(C.begin(), C.end(), v) != C.end());
      bool inTree = (find(tree.begin(), tree.end(), v) != tree.end());

      if (!(inClique && !inTree)) {
        next.push_back(v);
      }
    }

    filtered.swap(next);
  }

  return filtered;
}

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
      foundCliques.push_back(R);
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

  // P empty but X non-empty: R is NOT maximal, prune
  if (isEmpty(P))
    return;

  // Choose pivot from P ∪ X
  ui pivot = choosePivot(P, X);

  if (isPSubsetOfFoundClique(P)) {
    redendantChecks.push_back(P);
    redundancy++;
  }

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

  bronKerboschRecursive(R, P, X);

  cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
  cout << "Maximum Clique Size: " << maxCliqueSize << endl;

  cout << endl << endl << "Redundancy Count: " << redundancy << endl;

  cout << "Redundant Checks: ";
  for (const auto &check : redendantChecks) {
    cout << "{ ";
    for (ui v : check)
      cout << v << " ";
    cout << "} ";
  }
  cout << endl;
}

// Latest Reordering Bron-Kerbosch Implementation
ReorderBK2::ReorderBK2(Graph &g, bool sortMode) {
  graph = g;
  n = g.n;
  adjList.resize(n);
  adjList2.resize(n);
  cliqueCount = 0;
  maxCliqueSize = 0;

  if (sortMode) {
    order.resize(n);
    for (ui i = 0; i < n; i++)
      order[i] = i;

    sort(order.begin(), order.end(), [&](ui a, ui b) {
      if (graph.degree[a] != graph.degree[b])
        return graph.degree[a] > graph.degree[b]; // descending
      return a < b;                               // tie-break
    });

    // Maping between old IDs and new IDs
    old2new.resize(n);
    new2old.resize(n);

    for (ui newId = 0; newId < n; newId++) {
      ui oldId = order[newId];
      old2new[oldId] = newId;
      new2old[newId] = oldId;
    }

    // New CSR representation
    vector<ui> newDegree(n, 0);
    vector<ui> newOffset(n + 1, 0);
    vector<ui> newNeighbors;

    // count degrees
    for (ui oldU = 0; oldU < n; oldU++) {
      ui newU = old2new[oldU];
      newDegree[newU] = graph.degree[oldU];
    }

    // prefix sum
    for (ui i = 0; i < n; i++)
      newOffset[i + 1] = newOffset[i] + newDegree[i];

    newNeighbors.resize(newOffset[n]);

    vector<ui> cursor = newOffset;

    for (ui oldU = 0; oldU < n; oldU++) {
      ui newU = old2new[oldU];

      for (ui j = graph.offset[oldU]; j < graph.offset[oldU + 1]; j++) {
        ui oldV = graph.neighbors[j];
        ui newV = old2new[oldV];

        newNeighbors[cursor[newU]++] = newV;
      }
    }

    // replace graph with new CSR
    graph.degree = std::move(newDegree);
    graph.offset = std::move(newOffset);
    graph.neighbors = std::move(newNeighbors);

    cout << "Reordered Graph Adjacency Lists:" << endl;
    for (ui i = 0; i < n; i++) {
      for (ui j = graph.offset[old2new[i]]; j < graph.offset[old2new[i] + 1];
           j++) {
        cout << new2old[graph.neighbors[j]] << " ";
      }
    }
  }
  cout << endl;
  // Fill adjacency lists from graph
  for (ui i = 0; i < n; i++) {
    for (ui j = graph.offset[i]; j < graph.offset[i + 1]; j++) {
      ui neighbor = graph.neighbors[j];
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

  for (ui i = 0; i < n; i++) {
    cout << "Adjacency List for vertex " << new2old[i] << ": ";
    for (ui neighbor : adjList[i]) {
      cout << new2old[neighbor] << " ";
    }
    cout << endl;
  }

  cout << endl << "order neighbors only greater vertices:" << endl;
  for (ui i = 0; i < n; i++) {
    cout << "Adjacency List for vertex " << new2old[i] << ": ";
    for (ui neighbor : adjList2[i]) {
      cout << new2old[neighbor] << " ";
    }
    cout << endl;
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
      cout << new2old[v] << " ";
    }
    cout << "}" << endl;
  }
  cout << endl << "   ExpandTo: " << endl;
  for (vector<ui> tree : expandTo) {
    cout << "     { ";
    for (ui v : tree) {
      cout << new2old[v] << " ";
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
  // If yes, we do not need to continue further in this rCall, and trees will
  // be updated and new call will be made. If no, we continue to enemurate
  // trees in this rCall
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
    // Q = expandTo[i];
    Q = applyPreviousCliqueEffects(mustin[i], expandTo[i]);

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
    cout << new2old[v] << " ";
  }
  cout << endl;
  cout << "       Current Q: ";
  for (ui v : Q) {
    cout << new2old[v] << " ";
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

      // Set moveToNext to false to stop further rCall (tree expansions) as
      // tree need to be reordered
      moveToNext = false;
      if (debug) {
        cout << endl << "************************" << endl;
        cout << "Maximal Clique: { ";
        for (ui v : R) {
          cout << new2old[v] << " ";
        }
        cout << "}";
        cout << endl << "************************" << endl;
      }

      // Reordering logic to be added here

      // Create new trees (mustin and expandTo) for next rCall
      vector<vector<ui>> newMustin;
      vector<vector<ui>> newExpandTo;

      // removes trees who mustin vertices are all in maximal Clique(R).
      // For remaining trees, update their expandTo sets by including
      // Verticies in of maximal Clique (R), if they are connected to all
      // vertices in mustin.
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
          cout << new2old[v] << " ";
        }
        cout << "}" << endl;
      }
      cout << endl;

      cout << "       New ExpandTo after removing clique vertices: " << endl;
      for (vector<ui> tree : newExpandTo) {
        cout << "         { ";
        for (ui v : tree) {
          cout << new2old[v] << " ";
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
