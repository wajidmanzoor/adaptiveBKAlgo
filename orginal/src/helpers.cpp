#include "../inc/helpers.h"

void BronKerbosch(ull R, ull P, ull X) {
  byte i;
  ull v, single;
  if ((P == 0) && (X == 0)) {

    clique[count] = R;
    count++;
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
