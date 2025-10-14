#include "../inc/helpers.h"

void initializeAdjacencyMatrix(Graph& g) {
  // Initialize the adjacency matrix representation using bit masks
  N = g.n;
  if (N > Nmax) {
    std::cout << "Graph too large for bit representation (max " << Nmax << " vertices)" << std::endl;
    exit(1);
  }
  
  // Initialize pik array - pik[i] represents neighbors of vertex i as a bitmask
  for (int i = 0; i < N; i++) {
    pik[i] = 0;
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < N) {
        pik[i] |= (1ULL << neighbor);
      }
    }
  }
}

void BronKerbosch(ull R, ull P, ull X) {
  if ((P == 0) && (X == 0)) {
    // Found a maximal clique
    clique[::count] = R;
    
    // Print the clique
    std::cout << "Clique " << (::count + 1) << ": ";
    for (byte i = 0; i < N; i++) {
      if (R & (1ULL << i)) {
        std::cout << (int)i << " ";
      }
    }
    std::cout << std::endl;
    
    ::count++;
  } else {
    // Iterate through vertices in P
    for (byte i = 0; i < N; i++) {
      ull v = 1ULL << i;
      if (P & v) {
        // Recursive call: R ∪ {v}, P ∩ N(v), X ∩ N(v)
        BronKerbosch(R | v, P & pik[i], X & pik[i]);
        // Move v from P to X
        P = P ^ v;
        X = X | v;
      }
    }
  }
}

void BKstandard(Graph& g) {
  std::cout << "Starting Bron-Kerbosch algorithm..." << std::endl;
  std::cout << "Graph has " << g.n << " vertices and " << g.m << " edges" << std::endl;
  
  // Initialize global variables
  ::count = 0;
  
  // Initialize adjacency matrix representation
  initializeAdjacencyMatrix(g);
  
  // Initialize sets:
  // R = {} (empty set)
  // P = all vertices
  // X = {} (empty set)
  ull R = 0;
  ull P = (1ULL << N) - 1; // All vertices set to 1
  ull X = 0;
  
  
  // Start the algorithm
  BronKerbosch(R, P, X);
  
  std::cout << "Found " << (int)::count << " maximal cliques" << std::endl;
}

// Dynamic BK Implementation
DynamicBK::DynamicBK(Graph& g) : n(g.n), clique_count(0) {
    // Initialize adjacency matrix
    adjacency.resize(n, vector<bool>(n, false));
    
    // Fill adjacency matrix from graph
    for (ui i = 0; i < n; i++) {
        for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
            ui neighbor = g.neighbors[j];
            if (neighbor < n) {
                adjacency[i][neighbor] = true;
            }
        }
    }
    
    std::cout << "Dynamic BK initialized for " << n << " vertices" << std::endl;
}

vector<bool> DynamicBK::intersect(const vector<bool>& set1, const vector<bool>& set2) {
    vector<bool> result(n, false);
    for (ui i = 0; i < n; i++) {
        if (set1[i] && set2[i]) {
            result[i] = true;
        }
    }
    return result;
}

vector<bool> DynamicBK::setUnion(const vector<bool>& set1, const vector<bool>& set2) {
    vector<bool> result(n, false);
    for (ui i = 0; i < n; i++) {
        result[i] = set1[i] || set2[i];
    }
    return result;
}

vector<bool> DynamicBK::setDifference(const vector<bool>& set1, const vector<bool>& set2) {
    vector<bool> result(n, false);
    for (ui i = 0; i < n; i++) {
        result[i] = set1[i] && !set2[i];
    }
    return result;
}

bool DynamicBK::isEmpty(const vector<bool>& set) {
    for (ui i = 0; i < n; i++) {
        if (set[i]) return false;
    }
    return true;
}

void DynamicBK::printSet(const vector<bool>& set, const string& name) {
    std::cout << name << ": ";
    for (ui i = 0; i < n; i++) {
        if (set[i]) {
            std::cout << i << " ";
        }
    }
    std::cout << std::endl;
}

void DynamicBK::bronKerboschRecursive(vector<bool>& R, vector<bool>& P, vector<bool>& X) {
    if (isEmpty(P) && isEmpty(X)) {
        // Found a maximal clique
        clique_count++;
        std::cout << "Clique " << clique_count << ": ";
        for (ui i = 0; i < n; i++) {
            if (R[i]) {
                std::cout << i << " ";
            }
        }
        std::cout << std::endl;
    } else {
        // Create a copy of P to iterate over
        vector<bool> P_copy = P;
        
        for (ui v = 0; v < n; v++) {
            if (P_copy[v]) {
                // Create R ∪ {v}
                vector<bool> new_R = R;
                new_R[v] = true;
                
                // Create P ∩ N(v) and X ∩ N(v)
                vector<bool> new_P = intersect(P, adjacency[v]);
                vector<bool> new_X = intersect(X, adjacency[v]);
                
                // Recursive call
                bronKerboschRecursive(new_R, new_P, new_X);
                
                // Move v from P to X
                P[v] = false;
                X[v] = true;
            }
        }
    }
}

void DynamicBK::findAllMaximalCliques() {
    std::cout << "Starting Dynamic Bron-Kerbosch algorithm..." << std::endl;
    std::cout << "Graph has " << n << " vertices" << std::endl;
    
    // Initialize sets
    vector<bool> R(n, false);  // Current clique (empty)
    vector<bool> P(n, true);   // All vertices as candidates
    vector<bool> X(n, false);  // Excluded set (empty)
    
    clique_count = 0;
    bronKerboschRecursive(R, P, X);
    
    std::cout << "Found " << clique_count << " maximal cliques" << std::endl;
}

// Sparse Dynamic BK Implementation (more memory efficient for sparse graphs)
SparseDynamicBK::SparseDynamicBK(Graph& g) : n(g.n), clique_count(0) {
    // Initialize adjacency lists
    adjacency.resize(n);
    
    // Fill adjacency lists from graph
    for (ui i = 0; i < n; i++) {
        for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
            ui neighbor = g.neighbors[j];
            if (neighbor < n) {
                adjacency[i].push_back(neighbor);
            }
        }
        // Sort for efficient intersection operations
        sort(adjacency[i].begin(), adjacency[i].end());
    }
    
    std::cout << "Sparse Dynamic BK initialized for " << n << " vertices" << std::endl;
}

vector<ui> SparseDynamicBK::intersect(const vector<ui>& set1, const vector<ui>& neighbors) {
    vector<ui> result;
    result.reserve(min(set1.size(), neighbors.size()));
    
    // Use two-pointer technique for sorted vectors
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

bool SparseDynamicBK::isEmpty(const vector<ui>& set) {
    return set.empty();
}

void SparseDynamicBK::bronKerboschRecursive(vector<ui>& R, vector<ui>& P, vector<ui>& X) {
    if (isEmpty(P) && isEmpty(X)) {
        // Found a maximal clique
        clique_count++;
        std::cout << "Clique " << clique_count << ": ";
        for (ui vertex : R) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    } else {
        // Create a copy of P to iterate over
        vector<ui> P_copy = P;
        
        for (ui v : P_copy) {
            // Create R ∪ {v}
            vector<ui> new_R = R;
            new_R.push_back(v);
            
            // Create P ∩ N(v) and X ∩ N(v)
            vector<ui> new_P = intersect(P, adjacency[v]);
            vector<ui> new_X = intersect(X, adjacency[v]);
            
            // Recursive call
            bronKerboschRecursive(new_R, new_P, new_X);
            
            // Move v from P to X
            P.erase(find(P.begin(), P.end(), v));
            X.push_back(v);
        }
    }
}

void SparseDynamicBK::findAllMaximalCliques() {
    std::cout << "Starting Sparse Dynamic Bron-Kerbosch algorithm..." << std::endl;
    std::cout << "Graph has " << n << " vertices" << std::endl;
    
    // Initialize sets
    vector<ui> R;  // Current clique (empty)
    vector<ui> P;  // All vertices as candidates
    vector<ui> X;  // Excluded set (empty)
    
    // Fill P with all vertices
    for (ui i = 0; i < n; i++) {
        P.push_back(i);
    }
    
    clique_count = 0;
    bronKerboschRecursive(R, P, X);
    
    std::cout << "Found " << clique_count << " maximal cliques" << std::endl;
}
