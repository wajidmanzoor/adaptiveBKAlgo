#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"

#include <cerrno>
#include <cstring>
#include <stdexcept>

int main(int argc, const char *argv[]) {
  if (argc < 3 || argc > 6) {
    exit(1);
  }
  ull parsed = 0;
  if (!parseUnsigned(argv[2], 6, parsed) || (parsed != 0 && parsed != 1)) {
    cerr << "Invalid mode. Use 0 (PivotBK) or 1 "
            "(Pure reoder 128+PXR+ET+Budget).\n";
    return 1;
  }
  const ui mode = static_cast<ui>(parsed);
  ui order = 1;
  if (argc > 3) {
    if (!parseUnsigned(argv[3], 2, parsed)) {
      cerr << "Invalid order. Use 0, 1, or 2.\n";
      return 1;
    }
    order = static_cast<ui>(parsed);
  }
  if (mode == 0 && argc > 4) {
    cerr << "PivotBK accepts only <graph> 0 [order].\n";
    return 1;
  }
  ui minCliqueSize = 3;
  if (argc > 4) {
    if (!parseUnsigned(argv[4], UINT_MAX, parsed) || parsed == 0) {
      cerr << "Invalid minCliqueSize. Use an integer in 1.." << UINT_MAX
           << ".\n";
      return 1;
    }
    minCliqueSize = static_cast<ui>(parsed);
  }
  ull budget = 10000;
  if (const char *environmentBudget = getenv("PURE_HITSET_BUDGET")) {
    if (!parseUnsigned(environmentBudget, numeric_limits<ull>::max(), budget)) {
      cerr << "Invalid PURE_HITSET_BUDGET. Use a nonnegative integer.\n";
      return 1;
    }
  }
  if (argc > 5 && !parseUnsigned(argv[5], numeric_limits<ull>::max(), budget)) {
    cerr << "Invalid budget. Use a nonnegative integer.\n";
    return 1;
  }

  Graph graph(argv[1]);
  const DegOrder degreeOrder = static_cast<DegOrder>(order);

  if (mode == 0) {
    cout << "Running PivotBK...\n";
    PivotBK algorithm(graph, degreeOrder);
    algorithm.findAllMaximalCliques();
    return 0;
  }

  cout << "Running Pure Reorder 128+PXR+ET+Budget"
       << " (budget=" << budget << ", minCliqueSize=" << minCliqueSize
       << ")...\n";
  ReorderSib algorithm(graph, degreeOrder, minCliqueSize);
  algorithm.setSolverWorkBudget(budget);
  algorithm.findAllMaximalCliquesPure();

  if (environmentFlagIsOne("VLDB_VALIDATION") ||
      environmentFlagIsOne("VLDB_PRINT_CLIQUES"))
    printCanonicalCliques(algorithm.getCliques());
  return 1;
}
