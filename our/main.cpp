#include "inc/config.h"
#include "inc/graph.h"
#include "inc/helpers.h"

#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>

namespace {

struct CliOptions {
  string graphPath;
  bool budgetEnabled = true;
  ull budget = pure_config::kDefaultBudget;
  ui minCliqueSize = 3;
  bool printCliques = false;
};

void printUsage(const char *program) {
  cout << "Usage: " << program << " GRAPH [options]\n\n"
       << "Exact AdaptiveBK/ReorderSib maximal-clique enumeration.\n\n"
       << "Options:\n"
       << "  --budget N|unlimited  Compatibility-work budget per seed solver\n"
       << "                        call (default: 1000)\n"
       << "  --min-clique-size N   Output threshold (default: 3)\n"
       << "  --print-cliques       Print canonical original vertex IDs\n"
       << "  -h, --help            Show this help\n";
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

    Graph graph(options.graphPath);
    ReorderSib reorder(graph, options.minCliqueSize);
    if (options.budgetEnabled)
      reorder.setSolverWorkBudget(options.budget);
    else
      reorder.clearSolverWorkBudget();

    const auto start = chrono::steady_clock::now();
    reorder.findAllMaximalCliquesPure();
    const auto stop = chrono::steady_clock::now();
    const double runtimeMs =
        chrono::duration<double, milli>(stop - start).count();
    cout << fixed << setprecision(3)
         << "reorder.runtime_ms=" << runtimeMs << '\n';

    if (options.printCliques)
      printCanonicalCliques(reorder.getCliques());
    return 0;
  } catch (const std::exception &error) {
    cerr << "error: " << error.what() << '\n';
    return 2;
  }
}
