#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

int main(int argc, const char *argv[]) {
  if (argc != 2) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }

  string filepath = argv[1];

  Graph g(filepath);
  BKstandard(g);

  return 0;
}