#include "inc/early_termination.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hbbmc_faithful {

namespace {

using CliqueOptions = std::vector<std::vector<int>>;

std::uint64_t checked_add(std::uint64_t a, std::uint64_t b) {
  if (b > std::numeric_limits<std::uint64_t>::max() - a) {
    throw std::overflow_error("t-plex output count exceeds uint64_t");
  }
  return a + b;
}

std::uint64_t checked_multiply(std::uint64_t a, std::uint64_t b) {
  if (a != 0U && b > std::numeric_limits<std::uint64_t>::max() / a) {
    throw std::overflow_error("t-plex output count exceeds uint64_t");
  }
  return a * b;
}

// Number of ways to complete a maximal independent set of a path after the
// vertex at `last` has been selected. This is the count-only counterpart of
// path_rec and uses its F(r)=F(r-2)+F(r-3) recurrence in linear time.
std::uint64_t path_completion_count(std::size_t path_size, std::size_t last) {
  if (last >= path_size) {
    throw std::logic_error("selected path position is out of range");
  }
  const std::size_t remaining = path_size - last - 1U;
  std::vector<std::uint64_t> count(remaining + 1U, 1U);
  if (remaining >= 2U) {
    count[2] = 1U;
  }
  for (std::size_t r = 3U; r <= remaining; ++r) {
    count[r] = checked_add(count[r - 2U], count[r - 3U]);
  }
  return count[remaining];
}

std::uint64_t path_option_count(std::size_t path_size) {
  if (path_size == 0U) {
    return 1U;
  }
  std::uint64_t count = path_completion_count(path_size, 0U);
  if (path_size >= 2U) {
    count = checked_add(count, path_completion_count(path_size, 1U));
  }
  return count;
}

std::uint64_t cycle_option_count(std::size_t cycle_size) {
  if (cycle_size == 3U) {
    return 3U;
  }
  if (cycle_size == 4U) {
    return 2U;
  }
  if (cycle_size == 5U) {
    return 5U;
  }
  if (cycle_size < 3U) {
    throw std::logic_error("cycle component has fewer than three vertices");
  }
  const std::uint64_t first_two_cases =
      checked_multiply(2U, path_completion_count(cycle_size - 1U, 0U));
  const std::uint64_t third_case = path_completion_count(cycle_size - 4U, 0U);
  return checked_add(first_two_cases, third_case);
}

std::vector<std::vector<int>>
build_complement(const Graph &graph, const std::vector<int> &candidates) {
  const std::size_t n = candidates.size();
  std::vector<std::vector<int>> complement(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (!graph.adjacent(candidates[i], candidates[j])) {
        complement[i].push_back(static_cast<int>(j));
        complement[j].push_back(static_cast<int>(i));
      }
    }
  }
  return complement;
}

void path_rec(const std::vector<int> &path, std::size_t last,
              std::vector<int> &selected, CliqueOptions &output) {
  if (last + 2U >= path.size()) {
    // Copy element-by-element. Besides making the materialization point
    // explicit, this avoids a GCC 13 false-positive -Warray-bounds warning
    // caused by inlining vector's bulk-copy path through this recursion.
    output.emplace_back();
    output.back().reserve(selected.size());
    for (const int v : selected) {
      output.back().push_back(v);
    }
    return;
  }

  selected.push_back(path[last + 2U]);
  path_rec(path, last + 2U, selected, output);
  selected.pop_back();

  if (last + 3U < path.size()) {
    selected.push_back(path[last + 3U]);
    path_rec(path, last + 3U, selected, output);
    selected.pop_back();
  }
}

CliqueOptions path_options(const std::vector<int> &path) {
  if (path.empty()) {
    return {{}};
  }
  CliqueOptions output;
  std::vector<int> selected{path[0]};
  path_rec(path, 0U, selected, output);
  if (path.size() >= 2U) {
    selected.assign(1U, path[1]);
    path_rec(path, 1U, selected, output);
  }
  return output;
}

CliqueOptions path_options_with_initial(const std::vector<int> &path,
                                        std::vector<int> selected) {
  if (path.empty()) {
    return {std::move(selected)};
  }
  CliqueOptions output;
  path_rec(path, 0U, selected, output);
  return output;
}

CliqueOptions cycle_options(const std::vector<int> &cycle) {
  const std::size_t n = cycle.size();
  if (n == 3U) {
    return {{cycle[0]}, {cycle[1]}, {cycle[2]}};
  }
  if (n == 4U) {
    return {{cycle[0], cycle[2]}, {cycle[1], cycle[3]}};
  }
  if (n == 5U) {
    return {{cycle[0], cycle[2]},
            {cycle[0], cycle[3]},
            {cycle[1], cycle[3]},
            {cycle[1], cycle[4]},
            {cycle[2], cycle[4]}};
  }
  if (n < 3U) {
    throw std::logic_error("cycle component has fewer than three vertices");
  }

  CliqueOptions output;

  std::vector<int> case_one(cycle.begin(), cycle.end() - 1);
  auto choices = path_options_with_initial(case_one, {cycle[0]});
  output.insert(output.end(), choices.begin(), choices.end());

  std::vector<int> case_two(cycle.begin() + 1, cycle.end());
  choices = path_options_with_initial(case_two, {cycle[1]});
  output.insert(output.end(), choices.begin(), choices.end());

  std::vector<int> case_three(cycle.begin() + 2, cycle.end() - 2);
  choices = path_options_with_initial(case_three, {cycle.back(), cycle[2]});
  output.insert(output.end(), choices.begin(), choices.end());
  return output;
}

std::vector<int>
ordered_component(const std::vector<int> &component,
                  const std::vector<std::vector<int>> &complement,
                  const std::vector<int> &candidates, bool cycle) {
  int start = component.front();
  if (cycle) {
    start = *std::min_element(component.begin(), component.end(),
                              [&](int a, int b) {
                                return candidates[static_cast<std::size_t>(a)] <
                                       candidates[static_cast<std::size_t>(b)];
                              });
  } else {
    bool found = false;
    for (const int v : component) {
      if (complement[static_cast<std::size_t>(v)].size() == 1U &&
          (!found || candidates[static_cast<std::size_t>(v)] <
                         candidates[static_cast<std::size_t>(start)])) {
        start = v;
        found = true;
      }
    }
    if (!found) {
      throw std::logic_error("non-cycle complement component has no endpoint");
    }
  }

  std::vector<int> order;
  order.reserve(component.size());
  int previous = -1;
  int current = start;
  while (true) {
    order.push_back(candidates[static_cast<std::size_t>(current)]);
    int next = -1;
    for (const int neighbor : complement[static_cast<std::size_t>(current)]) {
      if (neighbor == previous) {
        continue;
      }
      if (cycle && previous == -1 && next != -1 &&
          candidates[static_cast<std::size_t>(neighbor)] >=
              candidates[static_cast<std::size_t>(next)]) {
        continue;
      }
      next = neighbor;
      if (!cycle || previous != -1) {
        break;
      }
    }
    if (next == -1 || (cycle && next == start)) {
      break;
    }
    previous = current;
    current = next;
    if (order.size() > component.size()) {
      throw std::logic_error("failed to order complement component");
    }
  }
  if (order.size() != component.size()) {
    throw std::logic_error("ordered complement component has wrong size");
  }
  return order;
}

void enumerate_product(const std::vector<CliqueOptions> &components,
                       std::size_t index, std::vector<int> &selected,
                       const CandidateCliqueCallback &callback) {
  if (index == components.size()) {
    callback(selected);
    return;
  }
  for (const auto &option : components[index]) {
    const std::size_t old_size = selected.size();
    selected.insert(selected.end(), option.begin(), option.end());
    enumerate_product(components, index + 1U, selected, callback);
    selected.resize(old_size);
  }
}

} // namespace

bool is_t_plex(const Graph &graph, const std::vector<int> &candidates, int t) {
  if (t < 1 || t > 3) {
    throw std::invalid_argument("t-plex threshold must be in [1,3]");
  }
  const int allowed_non_neighbors = t - 1;
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    int missing = 0;
    for (std::size_t j = 0; j < candidates.size(); ++j) {
      if (i != j && !graph.adjacent(candidates[i], candidates[j])) {
        ++missing;
        if (missing > allowed_non_neighbors) {
          return false;
        }
      }
    }
  }
  return true;
}

std::uint64_t count_t_plex_maximal_cliques(const Graph &graph,
                                           const std::vector<int> &candidates,
                                           int t) {
  if (!is_t_plex(graph, candidates, t)) {
    throw std::invalid_argument("candidate set is not the requested t-plex");
  }
  if (candidates.empty() || t == 1) {
    return 1U;
  }

  const auto complement = build_complement(graph, candidates);
  std::vector<unsigned char> visited(candidates.size(), 0U);
  std::uint64_t total = 1U;
  for (std::size_t source = 0; source < candidates.size(); ++source) {
    if (visited[source] != 0U) {
      continue;
    }
    if (complement[source].empty()) {
      visited[source] = 1U;
      continue;
    }

    std::vector<int> stack{static_cast<int>(source)};
    std::vector<int> component;
    visited[source] = 1U;
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      component.push_back(v);
      for (const int neighbor : complement[static_cast<std::size_t>(v)]) {
        if (visited[static_cast<std::size_t>(neighbor)] == 0U) {
          visited[static_cast<std::size_t>(neighbor)] = 1U;
          stack.push_back(neighbor);
        }
      }
    }

    std::uint64_t component_count = 0U;
    if (t == 2) {
      if (component.size() != 2U) {
        throw std::logic_error("2-plex complement is not a matching");
      }
      component_count = 2U;
    } else {
      bool cycle = true;
      for (const int v : component) {
        const std::size_t degree =
            complement[static_cast<std::size_t>(v)].size();
        if (degree > 2U) {
          throw std::logic_error("3-plex complement degree exceeds two");
        }
        if (degree != 2U) {
          cycle = false;
        }
      }
      component_count = cycle ? cycle_option_count(component.size())
                              : path_option_count(component.size());
    }
    total = checked_multiply(total, component_count);
  }
  return total;
}

void enumerate_t_plex_maximal_cliques(const Graph &graph,
                                      const std::vector<int> &candidates, int t,
                                      const CandidateCliqueCallback &callback) {
  if (!is_t_plex(graph, candidates, t)) {
    throw std::invalid_argument("candidate set is not the requested t-plex");
  }
  if (candidates.empty()) {
    callback({});
    return;
  }
  if (t == 1) {
    callback(candidates);
    return;
  }

  const auto complement = build_complement(graph, candidates);
  std::vector<unsigned char> visited(candidates.size(), 0U);
  std::vector<int> fixed;
  std::vector<CliqueOptions> component_options;

  for (std::size_t source = 0; source < candidates.size(); ++source) {
    if (visited[source] != 0U) {
      continue;
    }
    if (complement[source].empty()) {
      visited[source] = 1U;
      fixed.push_back(candidates[source]);
      continue;
    }

    std::vector<int> stack{static_cast<int>(source)};
    std::vector<int> component;
    visited[source] = 1U;
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      component.push_back(v);
      for (const int neighbor : complement[static_cast<std::size_t>(v)]) {
        if (visited[static_cast<std::size_t>(neighbor)] == 0U) {
          visited[static_cast<std::size_t>(neighbor)] = 1U;
          stack.push_back(neighbor);
        }
      }
    }

    if (t == 2) {
      if (component.size() != 2U) {
        throw std::logic_error("2-plex complement is not a matching");
      }
      component_options.push_back(
          {{candidates[static_cast<std::size_t>(component[0])]},
           {candidates[static_cast<std::size_t>(component[1])]}});
      continue;
    }

    bool is_cycle = true;
    for (const int v : component) {
      const std::size_t degree = complement[static_cast<std::size_t>(v)].size();
      if (degree > 2U) {
        throw std::logic_error("3-plex complement degree exceeds two");
      }
      if (degree != 2U) {
        is_cycle = false;
      }
    }
    const auto order =
        ordered_component(component, complement, candidates, is_cycle);
    component_options.push_back(is_cycle ? cycle_options(order)
                                         : path_options(order));
  }

  enumerate_product(component_options, 0U, fixed, callback);
}

} // namespace hbbmc_faithful
