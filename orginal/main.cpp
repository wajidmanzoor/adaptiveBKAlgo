#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc != 3) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }

  string filepath = argv[1];
  int mode = atoi(argv[2]);

  Graph g(filepath);

  if (mode == 0) {
    cout << "Running Pivot BK (Original order)..." << endl;
    PivotBK pivotBk(g, DegOrder::ORIGINAL);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 1) {
    cout << "Running Pivot BK (Ascending degeneracy)..." << endl;
    PivotBK pivotBk(g, DegOrder::ASCENDING);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 2) {
    cout << "Running Pivot BK (Descending degeneracy)..." << endl;
    PivotBK pivotBk(g, DegOrder::DESCENDING);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 3) {
    cout << "Running Reorder Bron-Kerbosch Algorithm..." << endl;
    Reorder reorder(g, DegOrder::ORIGINAL);
    reorder.findAllMaximalCliques();
  } else if (mode == 4) {
    cout << "Running Reorder (Ascending Degeneracy) Bron-Kerbosch "
            "Algorithm..."
         << endl;
    Reorder reorder(g, DegOrder::ASCENDING);
    reorder.findAllMaximalCliques();
  } else if (mode == 5) {
    cout << "Running Reorder (Descending Degeneracy) Bron-Kerbosch "
            "Algorithm..."
         << endl;
    Reorder reorder(g, DegOrder::DESCENDING);
    reorder.findAllMaximalCliques();
  } else if (mode == 6) {
    cout << "Running ReorderSib (Brute Force By Size) Algorithm..." << endl;
    ReorderSib reorder(g, DegOrder::ORIGINAL, SibMethod::BRUTE_FORCE);
    reorder.findAllMaximalCliques();
  } else if (mode == 7) {
    cout << "Running ReorderSib (Backtracking Branch And Bound) Algorithm..."
         << endl;
    ReorderSib reorder(g, DegOrder::ORIGINAL, SibMethod::BACKTRACKING);
    reorder.findAllMaximalCliques();
  } else if (mode == 8) {
    cout << "Running ReorderSib (Greedy Approximation) Algorithm..." << endl;
    ReorderSib reorder(g, DegOrder::ORIGINAL, SibMethod::GREEDY);
    reorder.findAllMaximalCliques();
  } else if (mode == 9) {
    cout << "Running ReorderSib (Bitmask Exact Search) Algorithm..." << endl;
    ReorderSib reorder(g, DegOrder::ORIGINAL, SibMethod::BITMASK);
    reorder.findAllMaximalCliques();
  } else if (mode == 10) {
    cout << "Running ReorderSib (Inclusion-Minimal Clique Hitting Set) "
            "Algorithm..."
         << endl;
    ReorderSib reorder(g, DegOrder::ORIGINAL, SibMethod::MIN_HITTING_SET);
    reorder.findAllMaximalCliques();
  }

  else {
    cout << "Invalid mode! Use 0..10." << endl;
    exit(1);
  }

  return 0;
}
