#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <graph_file>" << endl;
    exit(1);
  }

  string filepath = argv[1];
  
  // Load the graph
  Graph g(filepath);
  
  // Run the classical Bron-Kerbosch algorithm
  BKstandard(g);

  return 0;
}