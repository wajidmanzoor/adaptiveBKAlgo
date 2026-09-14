#include "../inc/graph.h"
#include "../inc/helpers.h"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <streambuf>
#include <utility>

namespace {

using CliqueSet = std::set<std::vector<ui>>;
constexpr ui kMaximumTestVertices = 12;

class NullBuffer final : public std::streambuf {
protected:
  int overflow(int character) override { return character; }
};

CliqueSet bruteForce(
    ui n, const bool edges[kMaximumTestVertices][kMaximumTestVertices],
    ui minimumSize) {
  CliqueSet result;
  const std::uint64_t subsetCount = std::uint64_t{1} << n;
  for (std::uint64_t mask = 0; mask < subsetCount; ++mask) {
    if (static_cast<ui>(__builtin_popcountll(mask)) < minimumSize)
      continue;

    bool clique = true;
    for (ui left = 0; left < n && clique; ++left) {
      if ((mask & (std::uint64_t{1} << left)) == 0)
        continue;
      for (ui right = left + 1; right < n; ++right) {
        if ((mask & (std::uint64_t{1} << right)) != 0 &&
            !edges[left][right]) {
          clique = false;
          break;
        }
      }
    }
    if (!clique)
      continue;

    bool maximal = true;
    for (ui outside = 0; outside < n && maximal; ++outside) {
      if ((mask & (std::uint64_t{1} << outside)) != 0)
        continue;
      bool extends = true;
      for (ui inside = 0; inside < n; ++inside) {
        if ((mask & (std::uint64_t{1} << inside)) != 0 &&
            !edges[outside][inside]) {
          extends = false;
          break;
        }
      }
      if (extends)
        maximal = false;
    }
    if (!maximal)
      continue;

    std::vector<ui> vertices;
    for (ui vertex = 0; vertex < n; ++vertex)
      if ((mask & (std::uint64_t{1} << vertex)) != 0)
        vertices.push_back(vertex);
    result.insert(std::move(vertices));
  }
  return result;
}

bool writeGraph(
    const char *path, ui n,
    const bool edges[kMaximumTestVertices][kMaximumTestVertices],
    ui edgeCount) {
  std::ofstream output(path, std::ios::out | std::ios::trunc);
  if (!output)
    return false;
  output << n << ' ' << edgeCount << '\n';
  for (ui vertex = 0; vertex < n; ++vertex) {
    output << vertex;
    for (ui neighbor = 0; neighbor < n; ++neighbor)
      if (edges[vertex][neighbor])
        output << ' ' << neighbor;
    output << '\n';
  }
  return static_cast<bool>(output);
}

bool checkConfiguration(
    const char *path, ui n,
    const bool edges[kMaximumTestVertices][kMaximumTestVertices],
    ui minimumSize, ull budget, std::uint64_t graphId) {
  Graph graph(path);
  ReorderSib algorithm(graph, minimumSize);
  algorithm.setSolverWorkBudget(budget);
  algorithm.findAllMaximalCliquesPure();

  const std::vector<std::vector<ui>> emitted = algorithm.getCliques();
  CliqueSet actual(emitted.begin(), emitted.end());
  const CliqueSet expected = bruteForce(n, edges, minimumSize);
  if (actual == expected && actual.size() == emitted.size())
    return true;

  std::fprintf(stderr,
               "FAIL n=%u graph_id=%llu min_size=%u budget=%llu "
               "expected=%zu emitted=%zu unique=%zu\n",
               n, static_cast<unsigned long long>(graphId), minimumSize,
               static_cast<unsigned long long>(budget), expected.size(),
               emitted.size(), actual.size());
  return false;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: ccr_exhaustive TEMP_GRAPH_PATH\n");
    return 2;
  }

  NullBuffer nullBuffer;
  std::streambuf *savedOutput = std::cout.rdbuf(&nullBuffer);
  std::uint64_t graphCount = 0;

  for (ui n = 0; n <= 6; ++n) {
    std::vector<std::pair<ui, ui>> possibleEdges;
    for (ui left = 0; left < n; ++left)
      for (ui right = left + 1; right < n; ++right)
        possibleEdges.emplace_back(left, right);

    const std::uint64_t graphMasks =
        std::uint64_t{1} << possibleEdges.size();
    for (std::uint64_t graphMask = 0; graphMask < graphMasks; ++graphMask) {
      bool edges[kMaximumTestVertices][kMaximumTestVertices] = {};
      ui edgeCount = 0;
      for (size_t edge = 0; edge < possibleEdges.size(); ++edge) {
        if ((graphMask & (std::uint64_t{1} << edge)) == 0)
          continue;
        const ui left = possibleEdges[edge].first;
        const ui right = possibleEdges[edge].second;
        edges[left][right] = edges[right][left] = true;
        ++edgeCount;
      }

      if (!writeGraph(argv[1], n, edges, edgeCount)) {
        std::cout.rdbuf(savedOutput);
        std::fprintf(stderr, "cannot write temporary graph: %s\n", argv[1]);
        return 2;
      }

      // Budget zero forces the full-CCRMCE fallback whenever a seed-solver call
      // performs compatibility work. The second configuration exercises the
      // normal budget and minimum-size pruning.
      if (!checkConfiguration(argv[1], n, edges, 1, 0, graphMask) ||
          !checkConfiguration(argv[1], n, edges, 3, 1000, graphMask)) {
        std::cout.rdbuf(savedOutput);
        return 1;
      }
      ++graphCount;
    }
  }

  // Exhausting every graph beyond six vertices is impractical. These fixed-
  // seed samples extend oracle comparison through larger Q states and span
  // very sparse, intermediate, and nearly complete induced graphs.
  std::mt19937_64 random(0x4343524d43455f31ULL);
  constexpr unsigned densitiesPerThousand[] = {50, 200, 500, 800, 950};
  std::uint64_t randomGraphCount = 0;
  for (ui n = 7; n <= kMaximumTestVertices; ++n) {
    for (unsigned density : densitiesPerThousand) {
      for (unsigned sample = 0; sample < 20; ++sample) {
        bool edges[kMaximumTestVertices][kMaximumTestVertices] = {};
        ui edgeCount = 0;
        for (ui left = 0; left < n; ++left) {
          for (ui right = left + 1; right < n; ++right) {
            if (random() % 1000 >= density)
              continue;
            edges[left][right] = edges[right][left] = true;
            ++edgeCount;
          }
        }

        if (!writeGraph(argv[1], n, edges, edgeCount)) {
          std::cout.rdbuf(savedOutput);
          std::fprintf(stderr, "cannot write temporary graph: %s\n", argv[1]);
          return 2;
        }

        const std::uint64_t graphId =
            (std::uint64_t{1} << 63) | randomGraphCount;
        if (!checkConfiguration(argv[1], n, edges, 1, 0, graphId) ||
            !checkConfiguration(argv[1], n, edges, 3, 1000, graphId)) {
          std::cout.rdbuf(savedOutput);
          return 1;
        }
        ++randomGraphCount;
        ++graphCount;
      }
    }
  }

  std::cout.rdbuf(savedOutput);
  std::remove(argv[1]);
  std::fprintf(stderr,
               "CCR_EXHAUSTIVE_PASS graphs=%llu configurations=%llu "
               "random_graphs=%llu max_vertices=%u\n",
               static_cast<unsigned long long>(graphCount),
               static_cast<unsigned long long>(graphCount * 2),
               static_cast<unsigned long long>(randomGraphCount),
               kMaximumTestVertices);
  return 0;
}
