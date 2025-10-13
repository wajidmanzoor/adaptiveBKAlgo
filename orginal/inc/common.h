#pragma once

#define miv(a, b) ((a) > (b) ? (b) : (a))
#define mav(a, b) ((a) < (b) ? (b) : (a))

#include <assert.h>
#include <string.h>

#include <cstdlib>
#include <fstream>
#include <list>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <limits.h>
#include <map>
#include <mutex>
#include <set>
#include <sys/stat.h>
#include <utility>
using namespace std;

typedef unsigned int ui;
typedef unsigned short ushort;
typedef unsigned char uchar;
#define Nmax 64 // Set the grid size limit

typedef unsigned char byte;
typedef unsigned long long ull;

const int INF = 1000000000;
const double DINF = 1e9;

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
