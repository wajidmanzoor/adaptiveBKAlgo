#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"
#include "inc/rmce_reduction.h"

#include <cerrno>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <stdexcept>

namespace {

bool environmentFlagIsOne(const char *name) {
  const char *value = std::getenv(name);
  return value != nullptr && std::strcmp(value, "1") == 0;
}

void printCanonicalClique(const vector<ui> &clique) {
  cout << "clique";
  for (ui vertex : clique)
    cout << ' ' << vertex;
  cout << '\n';
}

void printStoredCanonicalCliques(const vector<vector<ui>> &cliques) {
  for (vector<ui> clique : cliques) {
    sort(clique.begin(), clique.end());
    printCanonicalClique(clique);
  }
}

void printReducedCanonicalCliques(const RmceReductionResult &reduced,
                                  const vector<vector<ui>> &residualCliques) {
  printStoredCanonicalCliques(reduced.directlyEmittedCliques);
  for (vector<ui> clique : residualCliques) {
    for (ui &vertex : clique)
      vertex = reduced.residualToOriginal.at(vertex);
    sort(clique.begin(), clique.end());
    printCanonicalClique(clique);
  }
}

} // namespace

int runMain(int argc, const char *argv[]) {
  if (argc < 4 || argc > 14) {
    cout << "Usage: bk_algorithm <file> <mode> <meth> [hitSetLimit] "
            "[prune1] [prune2] [sp1] [sp2] [sp3] [sp4] [sp5] [sp6] "
            "[minCliqueSize]"
         << endl;
    cout << "  mode: 0=PivotBK  1=PureReorderExact" << endl;
    cout << "  order: ascending degeneracy (fixed)" << endl;
    cout << "  meth: 0=Backtracking  1=Optimized  (mode 1 only)"
         << endl;
    cout << "  mode 0: bk_algorithm <file> 0 <ignored-meth> [minCliqueSize]"
         << endl;
    cout << "  minCliqueSize: output threshold for both modes (default 3; "
            "use 1 for conventional MCE)"
         << endl;
    exit(1);
  }

  string filepath = argv[1];
  int mode = atoi(argv[2]);
  int meth = atoi(argv[3]);
  ui hitSetLimit = UINT_MAX;
  bool prune1 = true;
  bool prune2 = true;
  bool sp1 = true;
  bool sp2 = true;
  bool sp3 = true;
  bool sp4 = true;
  bool sp5 = true;
  bool sp6 = true;
  ui minCliqueSize = 3;

  auto parsePositiveUi = [](const char *text, const char *name,
                            ui &result) -> bool {
    char *end = nullptr;
    errno = 0;
    const unsigned long long parsed = strtoull(text, &end, 10);
    if (text[0] < '0' || text[0] > '9' || errno == ERANGE || end == text ||
        *end != '\0' || parsed == 0 || parsed > UINT_MAX) {
      cerr << "Invalid " << name << "! Use an integer in 1.." << UINT_MAX
           << "." << endl;
      return false;
    }
    result = static_cast<ui>(parsed);
    return true;
  };

  if (mode == 0) {
    if (argc > 5) {
      cerr << "Mode 0 accepts only: <file> 0 <ignored-meth> "
              "[minCliqueSize]"
           << endl;
      return 1;
    }
    if (argc == 5 &&
        !parsePositiveUi(argv[4], "minCliqueSize", minCliqueSize))
      return 1;
  } else if (mode == 1) {
    if (argc > 4 && !parsePositiveUi(argv[4], "hitSetLimit", hitSetLimit))
      return 1;
    prune1 = (argc > 5) ? (bool)atoi(argv[5]) : true;
    prune2 = (argc > 6) ? (bool)atoi(argv[6]) : true;
    sp1 = (argc > 7) ? (bool)atoi(argv[7]) : true;
    sp2 = (argc > 8) ? (bool)atoi(argv[8]) : true;
    sp3 = (argc > 9) ? (bool)atoi(argv[9]) : true;
    sp4 = (argc > 10) ? (bool)atoi(argv[10]) : true;
    sp5 = (argc > 11) ? (bool)atoi(argv[11]) : true;
    sp6 = (argc > 12) ? (bool)atoi(argv[12]) : true;
    if (argc > 13 &&
        !parsePositiveUi(argv[13], "minCliqueSize", minCliqueSize))
      return 1;
  } else {
    cerr << "Invalid mode! Use 0 or 1." << endl;
    return 1;
  }

  const bool printCliqueIdentities = environmentFlagIsOne("VLDB_VALIDATION") ||
                                     environmentFlagIsOne("VLDB_PRINT_CLIQUES");

  Graph g(filepath);

  if (mode == 0) {
    cout << "Running Pivot BK (Ascending degeneracy)..." << endl;
    PivotBK pivotBk(g, minCliqueSize);
    pivotBk.findAllMaximalCliques();
  } else if (mode == 1) {
    if (meth < 0 || meth > 1) {
      cout << "Invalid method! Use 0 (Backtracking) or 1 (Optimized)." << endl;
      exit(1);
    }

    {
      cout << "Running Pure ReorderSib (exact worklist) "
              "(Ascending degeneracy) ";

      if (meth == 0)
        cout << "(Backtracking Branch And Bound) Algorithm...";
      else
        cout << "(Optimized Exact Search) Algorithm...";
      cout << endl;

      RmceReductionResult reduced;
      Graph *searchGraph = &g;
      const bool useRmce = environmentFlagIsOne("PURE_RMCE");
      if (useRmce) {
        const auto reduceStart = chrono::steady_clock::now();
        reduced = applyRmceReduction(g, minCliqueSize, printCliqueIdentities);
        const auto reduceEnd = chrono::steady_clock::now();
        const double reductionMs =
            chrono::duration<double, milli>(reduceEnd - reduceStart).count();
        searchGraph = &reduced.graph;
        cout << fixed << setprecision(3)
             << "PureRMCE: direct=" << reduced.directlyEmittedCount
             << "  residualVertices=" << reduced.graph.n
             << "  residualEdges=" << reduced.graph.m
             << "  d0=" << reduced.counters.degree0Vertices
             << "  d1=" << reduced.counters.degree1Vertices
             << "  d2=" << reduced.counters.degree2Vertices
             << "  nontriangle=" << reduced.counters.nontriangleEdges
             << "  reductionTime=" << reductionMs << " ms" << endl;
      }

      ReorderSib reorder(*searchGraph, static_cast<SibMethod>(meth), hitSetLimit, prune1,
                         prune2, sp1, sp2, sp3, sp4, sp5, sp6, minCliqueSize);
      if (const char *budget = getenv("PURE_HITSET_BUDGET"))
        reorder.setSolverWorkBudget(strtoull(budget, nullptr, 10));
      if (useRmce)
        reorder.setExternalResults(reduced.directlyEmittedCount,
                                   reduced.maximumCliqueSize);
      reorder.findAllMaximalCliquesPure();
      if (printCliqueIdentities) {
        if (useRmce)
          printReducedCanonicalCliques(reduced, reorder.getCliques());
        else
          printStoredCanonicalCliques(reorder.getCliques());
      }
    }
  } else {
    cout << "Invalid mode! Use 0 or 1." << endl;
    exit(1);
  }

  return 0;
}

int main(int argc, const char *argv[]) {
  try {
    return runMain(argc, argv);
  } catch (const std::overflow_error &error) {
    cerr << "Counting overflow: " << error.what() << endl;
    return 2;
  }
}
