#include "../inc/common.h"

// Global variable definitions
byte N, count;
ull pik[Nmax], clique[Nmax], mask;

byte number(ull x) {
  byte i;
  ull single, ret;
  single = 1;
  ret = 0;
  for (i = 0; i < Nmax; i++) {
    ret = ret + ((x & single) >> i);
    single = single << 1;
  }
  return ret;
}