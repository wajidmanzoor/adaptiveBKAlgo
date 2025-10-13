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
        std::cout << i << " ";
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
