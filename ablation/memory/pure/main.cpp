#include "inc/config.h"
#include "inc/graph.h"
#include "inc/helpers.h"

#if !defined(PURE_LEAN_BENCHMARK)
#include <chrono>
#endif
#include <cerrno>
#include <cstdlib>
#include <stdexcept>

namespace {

struct CliOptions {
  string graphPath;
  SibMethod method = SibMethod::OPTIMIZED;
  bool budgetEnabled = true;
  ull budget = pure_config::kDefaultBudget;
  ui minCliqueSize = 3;
  ui smallQFullPxrThreshold = 4;
  bool printCliques = false;
};

void printUsage(const char *program) {
  cout << "Usage: " << program << " GRAPH [options]\n\n"
       << "Exact Pure ReorderSib maximal-clique enumeration.\n\n"
       << "Options:\n"
       << "  --method optimized|backtracking  Seed solver (default: optimized)\n"
       << "  --budget N|unlimited            Compatibility-work budget per seed\n"
       << "                                  solver call (default: 1000)\n"
       << "  --min-clique-size N             Output threshold (default: 3)\n"
       << "  --small-q-pxr N                Full PXR when |Q| <= N (0 disables)\n"
       << "                                  (experimental; default: 4)\n"
       << "  --print-cliques                 Print canonical original vertex IDs\n"
       << "  -h, --help                      Show this help\n";
}

bool parseUnsigned(const char *text, ull maximum, ull &value) {
  if (text == nullptr || text[0] < '0' || text[0] > '9')
    return false;
  char *end = nullptr;
  errno = 0;
  const unsigned long long parsed = strtoull(text, &end, 10);
  if (errno == ERANGE || end == text || *end != '\0' || parsed > maximum)
    return false;
  value = static_cast<ull>(parsed);
  return true;
}

bool nextValue(int argc, const char *argv[], int &index, const char *option,
               const char *&value) {
  if (++index >= argc) {
    cerr << "Missing value after " << option << ".\n";
    return false;
  }
  value = argv[index];
  return true;
}

bool parseCli(int argc, const char *argv[], CliOptions &options) {
  if (argc < 2) {
    printUsage(argv[0]);
    return false;
  }
  if (string(argv[1]) == "-h" || string(argv[1]) == "--help") {
    printUsage(argv[0]);
    exit(0);
  }

  options.graphPath = argv[1];
  for (int i = 2; i < argc; ++i) {
    const string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      exit(0);
    }
    if (arg == "--print-cliques") {
      options.printCliques = true;
      continue;
    }

    const char *value = nullptr;
    if (arg == "--method") {
      if (!nextValue(argc, argv, i, "--method", value))
        return false;
      const string method = value;
      if (method == "optimized")
        options.method = SibMethod::OPTIMIZED;
      else if (method == "backtracking")
        options.method = SibMethod::BACKTRACKING;
      else {
        cerr << "Invalid --method value: " << method << ".\n";
        return false;
      }
      continue;
    }
    if (arg == "--budget") {
      if (!nextValue(argc, argv, i, "--budget", value))
        return false;
      if (string(value) == "unlimited") {
        options.budgetEnabled = false;
        options.budget = 0;
      } else {
        ull parsed = 0;
        if (!parseUnsigned(value, numeric_limits<ull>::max(), parsed)) {
          cerr << "Invalid --budget value: " << value
               << ". Use a non-negative integer or unlimited.\n";
          return false;
        }
        options.budgetEnabled = true;
        options.budget = parsed;
      }
      continue;
    }
    if (arg == "--min-clique-size") {
      if (!nextValue(argc, argv, i, "--min-clique-size", value))
        return false;
      ull parsed = 0;
      if (!parseUnsigned(value, numeric_limits<ui>::max(), parsed) ||
          parsed == 0) {
        cerr << "Invalid --min-clique-size value: " << value << ".\n";
        return false;
      }
      options.minCliqueSize = static_cast<ui>(parsed);
      continue;
    }
    if (arg == "--small-q-pxr") {
      if (!nextValue(argc, argv, i, "--small-q-pxr", value))
        return false;
      ull parsed = 0;
      if (!parseUnsigned(value, numeric_limits<ui>::max(), parsed)) {
        cerr << "Invalid --small-q-pxr value: " << value << ".\n";
        return false;
      }
      options.smallQFullPxrThreshold = static_cast<ui>(parsed);
      continue;
    }


    cerr << "Unknown option: " << arg << ".\n";
    return false;
  }
  return true;
}

void printCanonicalCliques(const vector<vector<ui>> &cliques) {
  for (vector<ui> clique : cliques) {
    sort(clique.begin(), clique.end());
    cout << "clique";
    for (ui vertex : clique)
      cout << ' ' << vertex;
    cout << '\n';
  }
}

} // namespace

int main(int argc, const char *argv[]) {
  try {
    CliOptions options;
    if (!parseCli(argc, argv, options))
      return 2;
#if !defined(PURE_LEAN_BENCHMARK)
    const auto totalStart = chrono::high_resolution_clock::now();
#endif

    Graph graph(options.graphPath);
    ReorderSib pure(graph, options.method, options.minCliqueSize);
    pure.setSmallQFullPxrThreshold(options.smallQFullPxrThreshold);
    if (options.budgetEnabled)
      pure.setSolverWorkBudget(options.budget);
    else
      pure.clearSolverWorkBudget();

#if !defined(PURE_LEAN_BENCHMARK)
    const auto setupEnd = chrono::high_resolution_clock::now();
#endif
    pure.findAllMaximalCliquesPure();
#if !defined(PURE_LEAN_BENCHMARK)
    const auto totalEnd = chrono::high_resolution_clock::now();
    const double setupMs =
        chrono::duration<double, milli>(setupEnd - totalStart).count();
    const double totalMs =
        chrono::duration<double, milli>(totalEnd - totalStart).count();
    cout << "pure.setup_runtime_ms=" << setupMs << '\n'
         << "pure.total_runtime_ms=" << totalMs << '\n';
#endif
    if (options.printCliques)
      printCanonicalCliques(pure.getCliques());
    return 0;
  } catch (const std::exception &error) {
    cerr << "error: " << error.what() << '\n';
    return 2;
  }
}
