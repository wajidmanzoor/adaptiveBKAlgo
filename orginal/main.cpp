#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <graph_file> <mode>" << endl;
    cout << "Modes:" << endl;
    cout << "  0: Adjacency Matrix BK" << endl;
    cout << "  1: Adjacency List BK" << endl;
    cout << "  2: Pivot BK (Optimized)" << endl;
    cout << "  3: Reorder BK (Tree Reordering)" << endl;
    cout << "  4: Pure Reorder BK (Most Efficient)" << endl;
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
    cout << "Running Pivot Bron-Kerbosch (Optimized)..." << endl;
    PivotBK pivotBk(g);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 3) {
    cout << "Running Reorder Bron-Kerbosch (Tree Reordering)..." << endl;
    ReorderBK reorderBk(g);
    reorderBk.findAllMaximalCliques();
  } else if (mode == 4) {
    cout << "Running Pure Reorder Bron-Kerbosch (Most Efficient)..." << endl;
    PureReorderBK pureReorderBk(g);
    pureReorderBk.findAllMaximalCliques();
  } else {
    cout << "Invalid mode! Use 0, 1, 2, 3, or 4." << endl;
    exit(1);
  }

  return 0;
}