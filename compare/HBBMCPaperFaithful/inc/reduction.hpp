#pragma once

#include "graph.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace hbbmc_faithful {

struct ReductionOptions {
  // Count-only is the production default. Identity materialization and the
  // duplicate-detection set are enabled only for printing or validation.
  bool collect_cliques = false;
  bool validate_invariants = false;
  // Reporting filter only. RMCE reductions are applied identically.
  std::size_t minimum_output_clique_size = 1U;
};

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
// published HBBMC++.
class GraphReductionModule {
public:
  virtual ~GraphReductionModule() = default;
  virtual std::string name() const = 0;
  virtual bool complete_hbbmc_gr() const = 0;
  virtual ReductionResult apply(const Graph &graph,
                                ReductionOptions options = {}) const = 0;
};

std::unique_ptr<GraphReductionModule>
make_graph_reduction_module(const std::string &name);

} // namespace hbbmc_faithful
