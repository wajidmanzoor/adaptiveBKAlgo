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
#include <bitset>
using namespace std;

typedef unsigned int ui;
typedef unsigned short ushort;
typedef unsigned char uchar;
#define Nmax 64 // Set the grid size limit

typedef unsigned char byte;
typedef unsigned long long ull;

const int INF = 1000000000;
const double DINF = 1e9;

extern byte N, count;
extern ull pik[Nmax], clique[Nmax], mask;
byte number(ull x);
