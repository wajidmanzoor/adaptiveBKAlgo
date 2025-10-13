#include "../inc/helpers.h"

void intializeAdjacencyMatrix(Graph &g) {
  N = g.n;
  if (N > Nmax) {
    cout << "Graph too large for bit representation (max " << Nmax << " nodes)"
         << endl;
    exit(1);
  }

  for (int i = 0; i < N; i++) {
    pik[i] = 0;
    for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
      ui neighbor = g.neighbors[j];
      if (neighbor < N) {
        pik[i] |= (1ULL << neighbor);
      }
    }
  }

  for (int i = 0; i < N; i++) {
    cout << "pik[" << i << "]=" << bitset<Nmax>(pik[i]) << endl;
  }
}

void BronKerbosch(ull R, ull P, ull X) {
  cout << "R=" << bitset<Nmax>(R) << ", P=" << bitset<Nmax>(P)
       << ", X=" << bitset<Nmax>(X) << endl;
  if ((P == 0) && (X == 0)) {
    clique[::count++] = R;

    cout << "Clique " << (::count) << ": ";
    for (byte i = 0; i < N; i++) {
      if (R & (1ULL << i)) {
        cout << i << " ";
      }
    }
    cout << endl;
    ::count++;
  } else {
    for (byte i = 0; i < N; i++) {
      ull v = 1ULL << i;
      if (P & v) {
        BronKerbosch(R | v, P & pik[i], X & pik[i]);
        P = P ^ v;
        X = X | v;
      }
    }
  }
}

void BKstandard(Graph &g) {
  cout << "starting Standard Bron-Kerbosch algoritm" << endl;
  cout << "Graph has " << g.n << " nodes and " << g.m << " edges" << endl;

  ::count = 0;
  intializeAdjacencyMatrix(g);

  ull R = 0;
  ull P = (1ULL << N) - 1;
  ull X = 0;
  ::mask = (1ULL << N) - 1;

  BronKerbosch(R, P, X);
  cout << "Total cliques found: " << ::count << endl;
}