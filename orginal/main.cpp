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
    cout << "Running Pivot Bron-Kerbosch with pruning..." << endl;
    PivotBK pivotBk(g);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 3) {
    cout << "Running Depth-First Reordering Algorithm..." << endl;
    DepthFirstReorderBK depthFirstBk(g);
    depthFirstBk.findAllMaximalCliques();
  } else if (mode == 4) {
    cout << "Running New Reordering Bron-Kerbosch Algorithm..." << endl;
    ReorderBK reorderBk(g);
    reorderBk.findAllMaximalCliques();

  } else {
    cout << "Invalid mode! Use 0, 1, 2 or 3. " << endl;
    exit(1);
  }

  return 0;
}