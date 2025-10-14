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
    AdjMatBK adjMatBk(g);
    adjMatBk.findAllMaximalCliques();
  } else {
    AdjListBK adjListBk(g);
    adjListBk.findAllMaximalCliques();
  }

  return 0;
}