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
    cout << "Running Adjacency Matrix Bron-Kerbosch..." << endl;
    AdjMatBK adjMatBk(g);
    adjMatBk.findAllMaximalCliques();
  } else if (mode == 1) {
    cout << "Running Adjacency List Bron-Kerbosch..." << endl;
    AdjListBK adjListBk(g);
    adjListBk.findAllMaximalCliques();

  } else if (mode == 2) {
    cout << "Running Pivot BK (Original order)..." << endl;
    PivotBK pivotBk(g, DegOrder::ORIGINAL);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 3) {
    cout << "Running Pivot BK (Ascending degeneracy)..." << endl;
    PivotBK pivotBk(g, DegOrder::ASCENDING);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 8) {
    cout << "Running Pivot BK (Descending degeneracy)..." << endl;
    PivotBK pivotBk(g, DegOrder::DESCENDING);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 4) {
    cout << "Running Original Reordering Bron-Kerbosch Algorithm..." << endl;
    Reorder reorder(g);
    reorder.findAllMaximalCliques();
  } else if (mode == 5) {
    cout << "Running ReorderNew Bron-Kerbosch Algorithm..." << endl;
    ReorderNew reorderNew(g, DegOrder::ORIGINAL);
    reorderNew.findAllMaximalCliques();
  } else if (mode == 6) {
    cout << "Running ReorderNew (Ascending Degeneracy) Bron-Kerbosch "
            "Algorithm..."
         << endl;
    ReorderNew reorderNew(g, DegOrder::ASCENDING);
    reorderNew.findAllMaximalCliques();
  } else if (mode == 7) {
    cout << "Running ReorderNew (Descending Degeneracy) Bron-Kerbosch "
            "Algorithm..."
         << endl;
    ReorderNew reorderNew(g, DegOrder::DESCENDING);
    reorderNew.findAllMaximalCliques();
  }

  else {
    cout << "Invalid mode! Use 0, 1, 2, 5, 6, or 7." << endl;
    exit(1);
  }

  return 0;
}