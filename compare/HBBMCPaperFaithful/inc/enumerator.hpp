#pragma once

#include "graph.hpp"
#include "truss_order.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace hbbmc_faithful {

struct EnumerationOptions {
  // 0 gives bare HBBMC. Values 1..3 enable the paper's t-plex ET up to t.
  int early_termination_threshold = 0;
  // These two switches are jointly enabled by the complete "rmce" module.
  bool rmce_dynamic_reduction = false;
  bool rmce_forbidden_set_reduction = false;
  bool collect_cliques = false;
  bool validate_invariants = false;
  // Reporting filter only. Search, branching, and ET behavior are unchanged.
  std::size_t minimum_output_clique_size = 1U;
};

struct StructuralCounters {
  std::uint64_t root_edge_branches = 0;
  std::uint64_t descendant_edge_branches = 0;
  std::uint64_t root_candidate_vertices = 0;
  std::uint64_t root_excluded_vertices = 0;
  std::uint64_t root_candidate_max = 0;
  std::uint64_t root_excluded_max = 0;
  std::uint64_t root_support_mismatches = 0;
  std::uint64_t root_tau_bound_violations = 0;
  std::uint64_t root_partition_overlap_violations = 0;

  std::uint64_t vertex_recursive_calls = 0;
  std::uint64_t pivot_calls = 0;
  std::uint64_t pivot_branches = 0;
  std::uint64_t max_vertex_depth = 0;

  std::uint64_t dynamic_reduction_calls = 0;
  std::uint64_t dynamic_degree0_removed = 0;
  std::uint64_t dynamic_degree0_outputs = 0;
  std::uint64_t dynamic_degree1_removed = 0;
  std::uint64_t dynamic_degree1_outputs = 0;
  std::uint64_t dynamic_universal_moved = 0;
  std::uint64_t forbidden_set_removed = 0;

  std::uint64_t et_checks = 0;
  std::uint64_t et1_calls = 0;
  std::uint64_t et2_calls = 0;
  std::uint64_t et3_calls = 0;
  std::uint64_t et1_outputs = 0;
  std::uint64_t et2_outputs = 0;
  std::uint64_t et3_outputs = 0;

  std::uint64_t terminal_outputs = 0;
  std::uint64_t isolated_vertex_outputs = 0;
  std::uint64_t owner_violations = 0;
  std::uint64_t invalid_clique_outputs = 0;
  std::uint64_t nonmaximal_outputs = 0;
  std::uint64_t duplicate_outputs = 0;
};

struct EnumerationResult {
  std::uint64_t maximal_clique_count = 0;
  std::vector<std::vector<int>> cliques;
  TrussOrder truss_order;
  StructuralCounters counters;
};

class Enumerator {
public:
  Enumerator(const Graph &graph, EnumerationOptions options = {});
  EnumerationResult run();

private:
  void enumerate_root_edge(int edge_id, int rank);
  void vertex_bk(std::vector<int> &partial, std::vector<int> candidates,
                 std::vector<int> excluded, int owner_edge,
                 std::uint64_t depth);
  bool apply_rmce_dynamic_reduction(std::vector<int> &partial,
                                    std::vector<int> &candidates,
                                    std::vector<int> &excluded, int owner_edge);
  void apply_rmce_forbidden_set_reduction(const std::vector<int> &candidates,
                                          std::vector<int> &excluded);
  bool try_early_termination(const std::vector<int> &partial,
                             const std::vector<int> &candidates,
                             const std::vector<int> &excluded, int owner_edge);
  int choose_pivot(const std::vector<int> &candidates,
                   const std::vector<int> &excluded) const;
  bool edge_is_after_owner(int u, int v, int owner_edge) const;
  bool candidate_graph_matches_original(const std::vector<int> &candidates,
                                        int owner_edge) const;
  void record_output_count(std::size_t clique_size, std::uint64_t count = 1U);
  void emit_clique(std::vector<int> clique, int owner_edge);

  bool is_clique(const std::vector<int> &clique) const;
  bool is_maximal(const std::vector<int> &clique) const;
  std::string clique_key(const std::vector<int> &clique) const;

  const Graph &graph_;
  EnumerationOptions options_;
  EnumerationResult result_;
  std::unordered_set<std::string> seen_cliques_;
  std::vector<std::uint32_t> dynamic_mark_epoch_;
  std::uint32_t dynamic_mark_epoch_value_ = 0U;
};

} // namespace hbbmc_faithful
