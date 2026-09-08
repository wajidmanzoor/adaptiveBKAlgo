#pragma once

#include "graph.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace hbbmc_faithful {

struct ReductionOptions {
  bool validate_invariants = false;
  // Reporting filter only. RMCE reductions are applied identically.
  std::size_t minimum_output_clique_size = 1U;
};

// The outcome of a graph-reduction pass: a smaller ("residual") graph that
// the enumerator still needs to search, plus every maximal clique that the
// reduction rules could already prove complete on their own ("directly
// emitted"). directly_emitted_count + the enumerator's count on `graph`
// together equal the total number of retained maximal cliques in the input (see
// main.cpp's validate_complete_pipeline, and the README's "Together these
// facts make direct and residual outputs disjoint" argument).
struct ReductionResult {
  Graph graph;
  std::uint64_t directly_emitted_count = 0;
  std::vector<std::vector<std::int64_t>> directly_emitted_cliques;
  struct Counters {
    std::uint64_t degree0_vertices = 0;
    std::uint64_t degree1_vertices = 0;
    std::uint64_t degree2_vertices = 0;
    std::uint64_t nontriangle_edges = 0;
    std::uint64_t directly_emitted = 0;
    std::uint64_t duplicate_direct_outputs = 0;
  } counters;
};

// Deliberately modular boundary for the RMCE graph-reduction rules used by
// published HBBMC++. "none" is the identity module (bare HBBMC / HBBMC+ET);
// "rmce" (alias "gr") is the full global-reduction module implemented in
// reduction.cpp, whose apply() removes degree-<=2 vertices and
// triangle-free edges to a fixed point (RMCE Algorithms 3 and 4).
class GraphReductionModule {
public:
  virtual ~GraphReductionModule() = default;
  virtual std::string name() const = 0;
  // True only for the "rmce" module. main.cpp uses this to decide whether
  // to also turn on the enumerator's dynamic and forbidden-set reductions,
  // since those are only sound as a complete package (see the README's
  // "Feature and proof mapping" table).
  virtual bool complete_hbbmc_gr() const = 0;
  virtual ReductionResult apply(const Graph &graph,
                                ReductionOptions options = {}) const = 0;
};

std::unique_ptr<GraphReductionModule>
make_graph_reduction_module(const std::string &name);

} // namespace hbbmc_faithful
