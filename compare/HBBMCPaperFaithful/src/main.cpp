#include "inc/enumerator.hpp"
#include "inc/graph.hpp"
#include "inc/reduction.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct CliOptions {
  std::string graph_path;
  std::string graph_reduction = "none";
  int et_threshold = 0;
  int declared_vertices = -1;
  int minimum_clique_size = 1;
  bool print_cliques = false;
  bool print_counters = false;
  bool validate = false;
};

void usage(std::ostream &out, const char *program) {
  out << "Usage: " << program << " GRAPH [options]\n"
      << "\n"
      << "Paper-faithful clean-room reimplementation of HBBMC, HBBMC+,\n"
      << "and HBBMC++ (RMCE graph reduction plus t=3 early termination).\n"
      << "\nOptions:\n"
      << "  --et T                    early termination threshold 0..3 "
         "(default 0)\n"
      << "  --early-termination T     alias for --et\n"
      << "  --graph-reduction NAME   none or rmce (alias: gr)\n"
      << "  --min-clique-size N       report/print only cliques of size >= N\n"
      << "  --num-vertices N         vertex universe is 0..N-1 (retains "
         "isolates)\n"
      << "  --print-cliques          print each clique using original vertex "
         "labels\n"
      << "  --counters               print structural counters\n"
      << "  --validate               validate clique, maximality, uniqueness, "
         "ownership\n"
      << "  -h, --help               show this help\n";
}

int parse_int(const std::string &text, const std::string &option) {
  std::size_t used = 0;
  const int value = std::stoi(text, &used);
  if (used != text.size()) {
    throw std::invalid_argument("invalid integer for " + option + ": " + text);
  }
  return value;
}

CliOptions parse_cli(int argc, char **argv) {
  CliOptions options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage(std::cout, argv[0]);
      std::exit(0);
    }
    if (arg == "--print-cliques") {
      options.print_cliques = true;
      continue;
    }
    if (arg == "--counters") {
      options.print_counters = true;
      continue;
    }
    if (arg == "--validate") {
      options.validate = true;
      continue;
    }
    if (arg == "--et" || arg == "--early-termination" ||
        arg == "--graph-reduction" || arg == "--num-vertices" ||
        arg == "--min-clique-size") {
      if (++i >= argc) {
        throw std::invalid_argument("missing value after " + arg);
      }
      const std::string value = argv[i];
      if (arg == "--et" || arg == "--early-termination") {
        options.et_threshold = parse_int(value, arg);
      } else if (arg == "--graph-reduction") {
        options.graph_reduction = value;
      } else if (arg == "--min-clique-size") {
        options.minimum_clique_size = parse_int(value, arg);
      } else {
        options.declared_vertices = parse_int(value, arg);
      }
      continue;
    }
    if (!arg.empty() && arg[0] == '-') {
      throw std::invalid_argument("unknown option: " + arg);
    }
    if (!options.graph_path.empty()) {
      throw std::invalid_argument("more than one graph path supplied");
    }
    options.graph_path = arg;
  }
  if (options.graph_path.empty()) {
    throw std::invalid_argument("missing graph path");
  }
  if (options.et_threshold < 0 || options.et_threshold > 3) {
    throw std::invalid_argument("early-termination threshold must be in [0,3]");
    if (options.minimum_clique_size < 1) {
      throw std::invalid_argument("--min-clique-size must be positive");
    }
  }
  if (options.declared_vertices < -1) {
    throw std::invalid_argument("--num-vertices must be nonnegative");
  }
  return options;
}

void print_counter(const char *name, std::uint64_t value) {
  std::cout << "counter." << name << '=' << value << '\n';
}

void print_counters(const hbbmc_faithful::StructuralCounters &c) {
  print_counter("root_edge_branches", c.root_edge_branches);
  print_counter("descendant_edge_branches", c.descendant_edge_branches);
  print_counter("root_candidate_vertices", c.root_candidate_vertices);
  print_counter("root_excluded_vertices", c.root_excluded_vertices);
  print_counter("root_candidate_max", c.root_candidate_max);
  print_counter("root_excluded_max", c.root_excluded_max);
  print_counter("root_support_mismatches", c.root_support_mismatches);
  print_counter("root_tau_bound_violations", c.root_tau_bound_violations);
  print_counter("root_partition_overlap_violations",
                c.root_partition_overlap_violations);
  print_counter("vertex_recursive_calls", c.vertex_recursive_calls);
  print_counter("pivot_calls", c.pivot_calls);
  print_counter("pivot_branches", c.pivot_branches);
  print_counter("max_vertex_depth", c.max_vertex_depth);
  print_counter("dynamic_reduction_calls", c.dynamic_reduction_calls);
  print_counter("dynamic_degree0_removed", c.dynamic_degree0_removed);
  print_counter("dynamic_degree0_outputs", c.dynamic_degree0_outputs);
  print_counter("dynamic_degree1_removed", c.dynamic_degree1_removed);
  print_counter("dynamic_degree1_outputs", c.dynamic_degree1_outputs);
  print_counter("dynamic_universal_moved", c.dynamic_universal_moved);
  print_counter("forbidden_set_removed", c.forbidden_set_removed);
  print_counter("et_checks", c.et_checks);
  print_counter("et1_calls", c.et1_calls);
  print_counter("et2_calls", c.et2_calls);
  print_counter("et3_calls", c.et3_calls);
  print_counter("et1_outputs", c.et1_outputs);
  print_counter("et2_outputs", c.et2_outputs);
  print_counter("et3_outputs", c.et3_outputs);
  print_counter("terminal_outputs", c.terminal_outputs);
  print_counter("isolated_vertex_outputs", c.isolated_vertex_outputs);
  print_counter("owner_violations", c.owner_violations);
  print_counter("invalid_clique_outputs", c.invalid_clique_outputs);
  print_counter("nonmaximal_outputs", c.nonmaximal_outputs);
  print_counter("duplicate_outputs", c.duplicate_outputs);
}

void print_reduction_counters(
    const hbbmc_faithful::ReductionResult::Counters &c) {
  print_counter("gr_degree0_vertices", c.degree0_vertices);
  print_counter("gr_degree1_vertices", c.degree1_vertices);
  print_counter("gr_degree2_vertices", c.degree2_vertices);
  print_counter("gr_nontriangle_edges", c.nontriangle_edges);
  print_counter("gr_directly_emitted", c.directly_emitted);
  print_counter("gr_duplicate_direct_outputs", c.duplicate_direct_outputs);
}

void validate_complete_pipeline(
    const hbbmc_faithful::Graph &input,
    const hbbmc_faithful::ReductionResult &reduced,
    const hbbmc_faithful::EnumerationResult &residual) {
  std::unordered_map<std::int64_t, int> dense_by_label;
  dense_by_label.reserve(static_cast<std::size_t>(input.vertex_count()) * 2U +
                         1U);
  for (int v = 0; v < input.vertex_count(); ++v) {
    if (!dense_by_label.emplace(input.original_label(v), v).second) {
      throw std::logic_error("input contains duplicate original labels");
    }
  }

  std::set<std::vector<std::int64_t>> seen;
  auto validate_one = [&](std::vector<std::int64_t> labels) {
    std::sort(labels.begin(), labels.end());
    if (labels.empty() ||
        std::adjacent_find(labels.begin(), labels.end()) != labels.end()) {
      throw std::logic_error(
          "pipeline emitted an empty or repeated-label clique");
    }
    if (!seen.insert(labels).second) {
      throw std::logic_error("pipeline emitted a duplicate clique");
    }

    std::vector<int> clique;
    clique.reserve(labels.size());
    for (const std::int64_t label : labels) {
      const auto found = dense_by_label.find(label);
      if (found == dense_by_label.end()) {
        throw std::logic_error("pipeline emitted an unknown vertex label");
      }
      clique.push_back(found->second);
    }
    std::sort(clique.begin(), clique.end());
    for (std::size_t i = 0; i < clique.size(); ++i) {
      for (std::size_t j = i + 1U; j < clique.size(); ++j) {
        if (!input.adjacent(clique[i], clique[j])) {
          throw std::logic_error("pipeline emitted a non-clique");
        }
      }
    }
    for (int v = 0; v < input.vertex_count(); ++v) {
      if (std::binary_search(clique.begin(), clique.end(), v)) {
        continue;
      }
      bool extends = true;
      for (const int u : clique) {
        if (!input.adjacent(u, v)) {
          extends = false;
          break;
        }
      }
      if (extends) {
        throw std::logic_error("pipeline emitted a non-maximal clique");
      }
    }
  };

  for (const auto &clique : reduced.directly_emitted_cliques) {
    validate_one(clique);
  }
  for (const auto &clique : residual.cliques) {
    std::vector<std::int64_t> labels;
    labels.reserve(clique.size());
    for (const int v : clique) {
      labels.push_back(reduced.graph.original_label(v));
    }
    validate_one(std::move(labels));
  }
  if (seen.size() !=
      reduced.directly_emitted_count + residual.maximal_clique_count) {
    throw std::logic_error(
        "pipeline identity count differs from numeric count");
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    const CliOptions cli = parse_cli(argc, argv);
    const hbbmc_faithful::Graph input = hbbmc_faithful::Graph::read_edge_list(
        cli.graph_path, cli.declared_vertices);
    auto reduction =
        hbbmc_faithful::make_graph_reduction_module(cli.graph_reduction);

    const auto algorithm_start = std::chrono::steady_clock::now();
    const auto reduction_start = algorithm_start;
    hbbmc_faithful::ReductionOptions reduction_options;
    reduction_options.collect_cliques = cli.print_cliques || cli.validate;
    reduction_options.validate_invariants = cli.validate;
    reduction_options.minimum_output_clique_size =
        static_cast<std::size_t>(cli.minimum_clique_size);
    auto reduced = reduction->apply(input, reduction_options);
    const auto reduction_stop = std::chrono::steady_clock::now();

    hbbmc_faithful::EnumerationOptions options;
    options.early_termination_threshold = cli.et_threshold;
    options.rmce_dynamic_reduction = reduction->complete_hbbmc_gr();
    options.rmce_forbidden_set_reduction = reduction->complete_hbbmc_gr();
    options.collect_cliques = cli.print_cliques || cli.validate;
    options.validate_invariants = cli.validate;
    options.minimum_output_clique_size =
        static_cast<std::size_t>(cli.minimum_clique_size);

    const auto enumeration_start = std::chrono::steady_clock::now();
    hbbmc_faithful::Enumerator enumerator(reduced.graph, options);
    auto result = enumerator.run();
    const auto algorithm_stop = std::chrono::steady_clock::now();
    const double reduction_milliseconds =
        std::chrono::duration<double, std::milli>(reduction_stop -
                                                  reduction_start)
            .count();
    const double enumeration_milliseconds =
        std::chrono::duration<double, std::milli>(algorithm_stop -
                                                  enumeration_start)
            .count();
    const double total_milliseconds = std::chrono::duration<double, std::milli>(
                                          algorithm_stop - algorithm_start)
                                          .count();

    const std::uint64_t direct_count = reduced.directly_emitted_count;
    if (result.maximal_clique_count >
        std::numeric_limits<std::uint64_t>::max() - direct_count) {
      throw std::overflow_error("maximal-clique count exceeds uint64_t");
    }
    const std::uint64_t maximal_clique_count =
        direct_count + result.maximal_clique_count;
    if (cli.validate) {
      validate_complete_pipeline(input, reduced, result);
    }

    std::string algorithm;
    const bool complete_gr = reduction->complete_hbbmc_gr();
    const bool hbbmc_plus_plus =
        complete_gr && options.rmce_dynamic_reduction &&
        options.rmce_forbidden_set_reduction && cli.et_threshold == 3;
    if (hbbmc_plus_plus) {
      algorithm = "HBBMC++";
    } else if (complete_gr && cli.et_threshold == 0) {
      algorithm = "HBBMC+";
    } else if (complete_gr) {
      algorithm = "HBBMC+GR+ET(t=" + std::to_string(cli.et_threshold) + ")";
    } else if (cli.et_threshold == 0) {
      algorithm = "HBBMC";
    } else {
      algorithm = "HBBMC+ET(t=" + std::to_string(cli.et_threshold) + ")";
    }

    std::cout << "implementation=independent-paper-faithful-reimplementation\n";
    std::cout << "algorithm=" << algorithm << '\n';
    std::cout << "hbbmc_plus_plus=" << (hbbmc_plus_plus ? "true" : "false")
              << '\n';
    std::cout << "graph_reduction=" << reduction->name() << '\n';
    std::cout << "graph_reduction_complete=" << (complete_gr ? "true" : "false")
              << '\n';
    std::cout << "vertices=" << input.vertex_count() << '\n';
    std::cout << "edges=" << input.edge_count() << '\n';
    std::cout << "minimum_clique_size=" << cli.minimum_clique_size << '\n';
    std::cout << "residual_vertices=" << reduced.graph.vertex_count() << '\n';
    std::cout << "residual_edges=" << reduced.graph.edge_count() << '\n';
    std::cout << "truss_tau=" << result.truss_order.tau << '\n';
    std::cout << "direct_reduction_cliques=" << direct_count << '\n';
    std::cout << "residual_cliques=" << result.maximal_clique_count << '\n';
    std::cout << "maximal_cliques=" << maximal_clique_count << '\n';
    if (cli.validate) {
      std::cout << "pipeline_validation=pass\n";
    }
    std::cout << std::fixed << std::setprecision(3)
              << "reduction_runtime_ms=" << reduction_milliseconds << '\n'
              << "enumeration_runtime_ms=" << enumeration_milliseconds << '\n'
              << "algorithm_runtime_ms=" << total_milliseconds << '\n';

    if (cli.print_cliques) {
      for (const auto &clique : reduced.directly_emitted_cliques) {
        std::cout << "clique";
        for (const std::int64_t label : clique) {
          std::cout << ' ' << label;
        }
        std::cout << '\n';
      }
      for (const auto &clique : result.cliques) {
        std::cout << "clique";
        for (const int v : clique) {
          std::cout << ' ' << reduced.graph.original_label(v);
        }
        std::cout << '\n';
      }
    }
    if (cli.print_counters) {
      print_reduction_counters(reduced.counters);
      print_counters(result.counters);
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    usage(std::cerr, argv[0]);
    return 2;
  }
}
