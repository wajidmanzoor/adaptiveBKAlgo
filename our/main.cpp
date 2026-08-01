#include "inc/common.h"
#include "inc/fast_list_bk.h"
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
  if (argc < 5 || argc > 15) {
    cout << "Usage: bk_algorithm <file> <mode> <ord> <meth> [hitSetLimit] "
            "[prune1] [prune2] [sp1] [sp2] [sp3] [sp4] [sp5] [sp6] "
            "[minCliqueSize]"
         << endl;
    cout << "  mode: 0=PivotBK  1=HybridReorder  2=BitsetBK  3=LocalBitsetBK  "
            "4=AdaptiveBK  5=FastListBK  6=PureReorderExact"
         << endl;
    cout << "  ord:  0=Original  1=Ascending  2=Descending" << endl;
    cout << "  meth: 0=Backtracking  1=Optimized  (ReorderSib modes only)"
         << endl;
    cout << "  minCliqueSize: mode-1/5/6 output threshold (default 3; use 1 "
            "for conventional MCE)"
         << endl;
    exit(1);
  }

  string filepath = argv[1];
  int mode = atoi(argv[2]);
  int ord = atoi(argv[3]);
  int meth = atoi(argv[4]);
  ui hitSetLimit = UINT_MAX;
  if (argc > 5) {
    char *end = nullptr;
    errno = 0;
    const unsigned long long parsed = strtoull(argv[5], &end, 10);
    if (argv[5][0] < '0' || argv[5][0] > '9' || errno == ERANGE ||
        end == argv[5] || *end != '\0' || parsed == 0 || parsed > UINT_MAX) {
      cerr << "Invalid hitSetLimit! Use an integer in 1.." << UINT_MAX << "."
           << endl;
      return 1;
    }
    hitSetLimit = static_cast<ui>(parsed);
  }
  bool prune1 = (argc > 6) ? (bool)atoi(argv[6]) : true;
  bool prune2 = (argc > 7) ? (bool)atoi(argv[7]) : true;
  bool sp1 = (argc > 8) ? (bool)atoi(argv[8]) : true;
  bool sp2 = (argc > 9) ? (bool)atoi(argv[9]) : true;
  bool sp3 = (argc > 10) ? (bool)atoi(argv[10]) : true;
  bool sp4 = (argc > 11) ? (bool)atoi(argv[11]) : true;
  bool sp5 = (argc > 12) ? (bool)atoi(argv[12]) : true;
  bool sp6 = (argc > 13) ? (bool)atoi(argv[13]) : true;
  const bool minCliqueSizeExplicit = argc > 14;
  const bool printCliqueIdentities = environmentFlagIsOne("VLDB_VALIDATION") ||
                                     environmentFlagIsOne("VLDB_PRINT_CLIQUES");
  ui minCliqueSize = 3;
  if (argc > 14) {
    char *end = nullptr;
    errno = 0;
    const unsigned long long parsed = strtoull(argv[14], &end, 10);
    if (argv[14][0] < '0' || argv[14][0] > '9' || errno == ERANGE ||
        end == argv[14] || *end != '\0' || parsed == 0 || parsed > UINT_MAX) {
      cerr << "Invalid minCliqueSize! Use an integer in 1.." << UINT_MAX << "."
           << endl;
      return 1;
    }
    if (mode != 1 && mode != 5 && mode != 6) {
      cerr << "minCliqueSize is supported only by modes 1, 5, and 6; other "
              "modes retain their existing threshold."
           << endl;
      return 1;
    }
    minCliqueSize = static_cast<ui>(parsed);
  }

  Graph g(filepath);

  if (mode == 0) {
    cout << "Running Pivot BK ";
    if (ord == 0)
      cout << "(Original order)...";
    else if (ord == 1)
      cout << "(Ascending degeneracy)...";
    else if (ord == 2)
      cout << "(Descending degeneracy)...";
    else {
      cout << "Invalid order! Use 0..2." << endl;
      exit(1);
    }
    cout << endl;
    PivotBK pivotBk(g, static_cast<DegOrder>(ord));
    pivotBk.findAllMaximalCliques();
  } else if (mode == 2) {
    cout << "Running Bitset BK..." << endl;
    BitsetBK bitsetBk(g);
    bitsetBk.findAllMaximalCliques();
  } else if (mode == 3) {
    cout << "Running Local Bitset BK..." << endl;
    LocalBitsetBK localBitsetBk(g);
    localBitsetBk.findAllMaximalCliques();
  } else if (mode == 4) {
    const ull words = (g.n + 63) >> 6;
    const ull denseBytes = (ull)g.n * words * sizeof(ull);
    const double avgDegree =
        g.n == 0 ? 0.0 : (2.0 * static_cast<double>(g.m)) / g.n;
    ui maxDegree = 0;
    ui zeroDegree = 0;
    for (ui d : g.degree) {
      maxDegree = max(maxDegree, d);
      if (d == 0)
        zeroDegree++;
    }
    const bool cheapModuleDense = g.n <= 16000 && avgDegree >= 6.0 &&
                                  maxDegree >= 90 && zeroDegree * 4 >= g.n;
    const bool denseCoreCandidate = !cheapModuleDense &&
                                    denseBytes <= (64ULL << 20) &&
                                    avgDegree >= 3.0 && maxDegree >= 90;
    const ui degeneracy = denseCoreCandidate ? graphDegeneracy(g) : 0;
    const bool useDense =
        denseBytes <= (64ULL << 20) &&
        (g.n <= 1000 || (g.n <= 8000 && avgDegree >= 10.0) ||
         (g.n <= 10000 && avgDegree >= 14.0) || cheapModuleDense ||
         (g.n <= 16000 && avgDegree >= 6.0 && maxDegree >= 90 &&
          degeneracy >= 45) ||
         (g.n <= 16000 && avgDegree >= 40.0) ||
         (g.n <= 18000 && avgDegree >= 44.0) || avgDegree >= 46.0 ||
         degeneracy >= 64);

    if (useDense) {
      cout << "Running Adaptive BK (BitsetBK)..." << endl;
      BitsetBK bitsetBk(g);
      bitsetBk.findAllMaximalCliques();
    } else {
      cout << "Running Adaptive BK (LocalBitsetBK)..." << endl;
      LocalBitsetBK localBitsetBk(g);
      localBitsetBk.findAllMaximalCliques();
    }
  } else if (mode == 5) {
    cout << "Running Fast List BK..." << endl;
    FastListBK fastListBk(g, false, minCliqueSize);
    if (printCliqueIdentities)
      fastListBk.setCliqueSink(printCanonicalClique);
    fastListBk.findAllMaximalCliques();
  } else if (mode == 1 || mode == 6) {
    if (ord < 0 || ord > 2) {
      cout << "Invalid order! Use 0..2." << endl;
      exit(1);
    }
    if (meth < 0 || meth > 1) {
      cout << "Invalid method! Use 0 (Backtracking) or 1 (Optimized)." << endl;
      exit(1);
    }

    const bool allRulesEnabled =
        prune1 && prune2 && sp1 && sp2 && sp3 && sp4 && sp5 && sp6;
    const bool useAdaptivePivotExpansion =
        mode == 1 && ord == 1 && meth == 1 && hitSetLimit == UINT_MAX &&
        allRulesEnabled &&
        (g.n >= 128 || minCliqueSizeExplicit || printCliqueIdentities);

    if (mode == 1 && (minCliqueSizeExplicit || printCliqueIdentities) &&
        !useAdaptivePivotExpansion) {
      cerr << "An explicit minCliqueSize or clique-identity validation in "
              "mode 1 requires ord=1, meth=1, hitSetLimit=UINT_MAX, and all "
              "rule flags enabled so the Hybrid FastList lane is selected."
           << endl;
      return 1;
    }

    if (useAdaptivePivotExpansion) {
      cout << "Running ReorderSib (Hybrid reorder/sibling + pivot expansion)..."
           << endl;
      FastListBK fastListBk(g, true, minCliqueSize);
      if (printCliqueIdentities)
        fastListBk.setCliqueSink(printCanonicalClique);
      fastListBk.findAllMaximalCliques("ReorderSib");
    } else {
      // Mode 6 is the theorem-aligned Pure worklist.  Mode 1 retains the
      // original recursive implementation on configurations that do not route
      // to FastListBK, preserving the legacy lane for ablations.
      cout << (mode == 6 ? "Running Pure ReorderSib (exact worklist) "
                         : "Running ReorderSib Legacy ");
      if (ord == 0)
        cout << "(Original order) ";
      else if (ord == 1)
        cout << "(Ascending degeneracy) ";
      else
        cout << "(Descending degeneracy) ";

      if (meth == 0)
        cout << "(Backtracking Branch And Bound) Algorithm...";
      else
        cout << "(Optimized Exact Search) Algorithm...";
      cout << endl;

      RmceReductionResult reduced;
      Graph *searchGraph = &g;
      const bool useRmce = mode == 6 && environmentFlagIsOne("PURE_RMCE");
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

      ReorderSib reorder(*searchGraph, static_cast<DegOrder>(ord),
                         static_cast<SibMethod>(meth), hitSetLimit, prune1,
                         prune2, sp1, sp2, sp3, sp4, sp5, sp6, minCliqueSize);
      if (const char *budget = getenv("PURE_HITSET_BUDGET"))
        reorder.setSolverWorkBudget(strtoull(budget, nullptr, 10));
      if (useRmce)
        reorder.setExternalResults(reduced.directlyEmittedCount,
                                   reduced.maximumCliqueSize);
      if (mode == 6)
        reorder.findAllMaximalCliquesPure();
      else
        reorder.findAllMaximalCliques();
      if (mode == 6 && printCliqueIdentities) {
        if (useRmce)
          printReducedCanonicalCliques(reduced, reorder.getCliques());
        else
          printStoredCanonicalCliques(reorder.getCliques());
      }
    }
  } else {
    cout << "Invalid mode! Use 0..6." << endl;
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
