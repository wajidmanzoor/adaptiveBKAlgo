#include "../inc/common.h"
#include "../inc/graph.h"
#include "../inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <graph_file>" << endl;
    exit(1);
  }

  string filepath = argv[1];
  
  // Load graph
  Graph graph(filepath);
  
  // Convert to bitmap representation
  graph.convertToBitmap(pik, N);
  
  // Run classical BK algorithm
  BKstandard();
  
  // Print the found cliques
  printCliques();
  
  return 0;
}