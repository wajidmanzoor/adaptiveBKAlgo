#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc < 5 || argc > 14) {
    cout << "Usage: bk_algorithm <file> <mode> <ord> <meth> [hitSetLimit] [prune1] [prune2] [sp1] [sp2] [sp3] [sp4] [sp5] [sp6]" << endl;
    exit(1);
  }

  string filepath = argv[1];
  int mode = atoi(argv[2]);
  int ord  = atoi(argv[3]);
  int meth = atoi(argv[4]);
  ui   hitSetLimit = (argc > 5)  ? (ui)atoi(argv[5])    : UINT_MAX;
  bool prune1      = (argc > 6)  ? (bool)atoi(argv[6])  : true;
  bool prune2      = (argc > 7)  ? (bool)atoi(argv[7])  : true;
  bool sp1         = (argc > 8)  ? (bool)atoi(argv[8])  : true;
  bool sp2         = (argc > 9)  ? (bool)atoi(argv[9])  : true;
  bool sp3         = (argc > 10) ? (bool)atoi(argv[10]) : true;
  bool sp4         = (argc > 11) ? (bool)atoi(argv[11]) : true;
  bool sp5         = (argc > 12) ? (bool)atoi(argv[12]) : true;
  bool sp6         = (argc > 13) ? (bool)atoi(argv[13]) : true;

  Graph g(filepath);

  if (mode == 0) {
    cout << "Running Pivot BK ";
    if (ord == 0)
      cout << "(Original order)...";
    else if (ord == 1)
      cout << "(Ascending degeneracy)...";
    else if (ord == 2)
      cout << "(Descending degeneracy)...";
    else {
      cout << "Invalid order! Use 0..2." << endl;
      exit(1);
    }
    cout << endl;
    PivotBK pivotBk(g, static_cast<DegOrder>(ord));
    pivotBk.findAllMaximalCliques();
  } else if (mode == 1) {
    cout << "Running ReorderSib ";
    if (ord == 0)
      cout << "(Original order) ";
    else if (ord == 1)
      cout << "(Ascending degeneracy) ";
    else if (ord == 2)
      cout << "(Descending degeneracy) ";
    else {
      cout << "Invalid order! Use 0..2." << endl;
      exit(1);
    }
    if (meth == 0)
      cout << "(Brute Force By Size) Algorithm...";
    else if (meth == 1)
      cout << "(Backtracking Branch And Bound) Algorithm...";
    else if (meth == 2)
      cout << "(Greedy Approximation) Algorithm...";
    else if (meth == 3)
      cout << "(Bitmask Exact Search) Algorithm...";
    else if (meth == 4)
      cout << "(Inclusion-Minimal Clique Hitting Set) Algorithm...";
    else if (meth == 5)
      cout << "(Optimized Exact Search) Algorithm...";
    else {
      cout << "Invalid method! Use 0..5." << endl;
      exit(1);
    }
    cout << endl;
    ReorderSib reorder(g, static_cast<DegOrder>(ord),
                       static_cast<SibMethod>(meth), hitSetLimit,
                       prune1, prune2, sp1, sp2, sp3, sp4, sp5, sp6);
    reorder.findAllMaximalCliques();
  } else {
    cout << "Invalid mode! Use 0 or 1." << endl;
    exit(1);
  }

  return 0;
}
