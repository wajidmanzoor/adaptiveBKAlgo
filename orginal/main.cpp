#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc < 2 || argc > 3) {
    cout << "Usage: " << argv[0] << " <graph_file> [mode]" << endl;
    cout << "  mode: 'legacy' for bit-based (max 64 vertices, default)" << endl;
    cout << "        'dynamic' for vector-based (unlimited vertices)" << endl;
    cout << "        'sparse' for sparse adjacency lists (memory efficient)" << endl;
    exit(1);
  }

  string filepath = argv[1];
  string mode = "legacy";
  if (argc == 3) {
    mode = argv[2];
  }
  
  // Load the graph
  Graph g(filepath);
  
  if (mode == "dynamic") {
    // Use dynamic implementation (unlimited size)
    DynamicBK dynamicBK(g);
    dynamicBK.findAllMaximalCliques();
  } else if (mode == "sparse") {
    // Use sparse implementation (memory efficient)
    SparseDynamicBK sparseBK(g);
    sparseBK.findAllMaximalCliques();
  } else if (mode == "legacy") {
    // Use legacy bit-based implementation (max 64 vertices)
    if (g.n > 64) {
      cout << "Error: Legacy mode supports maximum 64 vertices. Graph has " << g.n << " vertices." << endl;
      cout << "Use 'dynamic' or 'sparse' mode for larger graphs." << endl;
      exit(1);
    }
    BKstandard(g);
  } else {
    cout << "Error: Unknown mode '" << mode << "'. Use 'legacy', 'dynamic', or 'sparse'." << endl;
    exit(1);
  }

  return 0;
}