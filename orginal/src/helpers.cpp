#include "../inc/helpers.h"

void BronKerbosch(ull R, ull P, ull X) {
  byte i;
  ull v, single;
  if ((P == 0) && (X == 0)) {
    clique[clique_count] = R;
    clique_count++;
  } else {
    single = 1;
    for (i = 0; i < N; i++) {
      v = single << (N - 1 - i);
      if (P & v) {
        P = P ^ v;
        BronKerbosch(R | v, P & pik[i], X & pik[i]);
        X = X | v;
      }
    }
  }
  return;
}

void BKstandard() {
  clique_count = 0;
  ull R = 0;  // Current clique
  ull P = (1ULL << N) - 1;  // All vertices initially in P
  ull X = 0;  // Empty X set
  
  BronKerbosch(R, P, X);
  
  cout << "Found " << (int)clique_count << " maximal cliques" << endl;
}

void printCliques() {
  cout << "Maximal cliques:" << endl;
  for (int i = 0; i < clique_count; i++) {
    cout << "Clique " << i + 1 << ": {";
    bool first = true;
    for (int j = 0; j < N; j++) {
      if (clique[i] & (1ULL << (N - 1 - j))) {
        if (!first) cout << ", ";
        cout << j;
        first = false;
      }
    }
    cout << "}" << endl;
  }
}
