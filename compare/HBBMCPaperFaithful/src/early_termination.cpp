#include "../inc/early_termination.hpp"

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

// Number of ways to complete a maximal independent set (MIS) of a path of
// `path_size` vertices, given that the vertex at index `last` was just
// selected (so index last+1 is already dominated/excluded). The remaining
// undecided suffix has length r = path_size-last-1. For r>=3 the same two
// choices as path_rec apply: select the next free vertex (leaves r-2
// undecided) or skip it and force-select the one after (leaves r-3
// undecided), giving count[r] = count[r-2] + count[r-3]. Base cases:
// count[0]=count[1]=1 (nothing left to decide, or one dominated leftover
// vertex -- exactly one way), count[2]=1 (two leftover vertices, only the
// "select the first of the two" branch is valid, matching path_rec's
// `last+3 < path.size()` guard). This is the count-only counterpart of
// path_rec, computed bottom-up in O(r) instead of by recursion.
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

// Total number of maximal independent sets of a path of path_size vertices.
// Mirrors path_options: the first selected vertex is either index 0 (then
// count via path_completion_count(.., 0)) or, when vertex 0 is left
// unselected, index 1 must be selected to dominate it (then count via
// path_completion_count(.., 1)); these two cases are disjoint and exhaustive
// for path_size>=2, so their counts simply add.
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

// Total number of maximal independent sets of a cycle of cycle_size
// vertices. Mirrors cycle_options's 3-way case split (see that function for
// the geometric argument): sizes 3/4/5 are hardcoded because the general
// recurrence below needs a path of at least 2 vertices for the "gap of two"
// case, and n=3 makes that path length negative/zero; for n>=6, "cycle[0]
// selected" and "cycle[1] selected instead" are symmetric and each reduce
// to a path completion count on an (n-1)-path (hence `2 *
// path_completion_count(cycle_size-1, 0)`), and the remaining "neither
// cycle[0] nor cycle[1] selected" case reduces to a path completion count
// on an (n-4)-path (the two forced-selected boundary vertices plus their
// two dominated neighbors are removed from consideration).
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

// Builds the complement of G[candidates] as an adjacency list indexed by
// *position in `candidates`* (not by graph vertex id): complement[i]
// contains j whenever candidates[i] and candidates[j] are NOT adjacent in
// the original graph. Since candidates is a t-plex, every vertex here has
// complement-degree <= t-1 <= 2, so this complement graph is a disjoint
// union of isolated vertices, paths, and (for t=3) cycles -- see the
// connected-component walks in count_t_plex_maximal_cliques and
// enumerate_t_plex_maximal_cliques below.
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

// Recursively enumerates every maximal independent set (MIS) of the path
// `path`, given that `path[last]` has already been chosen and appended to
// `selected` (so `path[last+1]`, its only still-relevant neighbor, is
// already dominated and must stay unselected). Each full assignment reached
// is appended to `output`. This is the direct-construction counterpart of
// path_completion_count's recurrence: from position `last`, either select
// `path[last+2]` (dominating `last+1` via `last`, and `last+2` itself
// selected) or, only when a vertex exists two further out, skip
// `path[last+2]` and force-select `path[last+3]` to dominate it instead.
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

// Every maximal independent set of `path` (a simple path given as an
// ordered vertex sequence), as full vertex-id lists. Same case split as
// path_option_count: either path[0] is selected (first path_rec call), or
// path[0] is left unselected and path[1] must be selected to dominate it
// (second call, only possible when the path has >=2 vertices).
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

// Like path_options, but `path[0]` is *not* itself part of the decision --
// the caller has already committed to some `selected` prefix (typically
// including path[0] as a forced choice) and only the rest of the path
// (positions 1 onward, since path_rec starts at last=0) is enumerated.
// Used by cycle_options to graft a path enumeration onto a cycle vertex
// that was fixed as selected ahead of time.
CliqueOptions path_options_with_initial(const std::vector<int> &path,
                                        std::vector<int> selected) {
  if (path.empty()) {
    return {std::move(selected)};
  }
  CliqueOptions output;
  path_rec(path, 0U, selected, output);
  return output;
}

// Every maximal independent set of `cycle` (a simple cycle given as an
// ordered vertex sequence, cycle[i] adjacent to cycle[i+1] and cycle[0]
// adjacent to cycle.back()), as full vertex-id lists. n=3/4/5 are
// hardcoded (small enough to enumerate by hand, and too small for the
// general reduction below to apply cleanly). For n>=6, every MIS falls
// into exactly one of three disjoint cases based on how cycle[0] gets
// dominated, each of which reduces to a path enumeration:
//   case_one:   cycle[0] is itself selected. Its two neighbors (cycle[1]
//               and cycle.back()) are then dominated "for free", so what
//               remains is exactly a path enumeration over cycle[0..n-2]
//               (cycle.back() dropped) with cycle[0] pre-selected.
//   case_two:   cycle[0] is not selected but cycle[1] is (dominating it).
//               By the mirror argument, remaining is a path enumeration
//               over cycle[1..n-1] with cycle[1] pre-selected.
//   case_three: neither cycle[0] nor cycle[1] is selected, so cycle[0]
//               must be dominated by cycle.back() and cycle[1] must be
//               dominated by cycle[2] -- both are therefore forced
//               selected, leaving a path enumeration over the untouched
//               middle segment cycle[2..n-3] with {cycle.back(), cycle[2]}
//               pre-selected.
// These three cases can't overlap (they differ on whether/which of
// cycle[0], cycle[1] is selected) and every MIS falls into one of them, so
// simple concatenation (not a cartesian product) gives the full result.
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

// Walks a connected component of the complement graph (given as unordered
// `component` positions, each with complement-degree 1 or 2) into an
// ordered sequence of *candidate vertex ids* following the path/cycle's
// actual edges, e.g. [c0, c1, c2, ...] with c_i adjacent to c_{i+1} in the
// complement graph (and, if `cycle`, c_last adjacent back to c0). This
// ordered form is what path_options/cycle_options above expect as input --
// they only know how to reason about "the vertex at position i" through
// its position in a linear/circular walk.
//
// Choice of starting vertex is deterministic (by original candidate vertex
// id) purely so the algorithm's output is reproducible across runs, not
// for correctness: cycles start at the globally minimum-labeled vertex
// (a cycle has no natural endpoint); paths start at whichever
// complement-degree-1 endpoint has the smaller label (a 2-vertex path has
// two, larger components have exactly two symmetric endpoints).
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

  // Walk forward: at each step, move to whichever complement-neighbor of
  // `current` isn't where we just came from (`previous`). A path/cycle
  // vertex has at most 2 complement-neighbors, so this is unambiguous
  // except at the very first step of a cycle walk (previous == -1, both
  // neighbors still "unvisited"), where the extra `candidates[neighbor] >=
  // candidates[next]` comparison picks the smaller-labeled of the two to
  // fix a single deterministic walk direction around the cycle.
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
    // A path ends when its far endpoint (degree 1, only neighbor already
    // visited) leaves no unvisited neighbor; a cycle ends when the walk
    // returns to `start`.
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

// Cartesian product over per-component MIS choices: since the complement
// graph's connected components share no edges, a full maximal independent
// set of the whole complement graph is exactly one MIS choice from each
// component combined together, independently. `selected` already holds the
// vertices forced by complement-isolated candidates (see `fixed` in
// enumerate_t_plex_maximal_cliques) before this recursion starts.
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

// True iff every candidate has at most t-1 non-neighbors among the other
// candidates, i.e. G[candidates] is a t-plex. Brute-force O(|candidates|^2)
// check; candidate sets reaching this function are the (small) contents of
// a BK recursion's P set, so this is cheap relative to continuing to branch.
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

// See the file-level count_t_plex_maximal_cliques doc comment in the header
// for its purpose (test-only algebraic cross-check). Structure mirrors
// enumerate_t_plex_maximal_cliques below but multiplies component sizes
// instead of building each combination, so it never touches path_rec /
// enumerate_product's O(output size) cost.
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
  // Depth-first walk (via explicit stack) over each connected component of
  // the complement graph; complement-isolated vertices (t-plex candidates
  // adjacent to every other candidate) contribute a factor of 1 and are
  // skipped, matching their forced-included treatment below.
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
      // A 2-plex's complement has max degree 1, so every non-trivial
      // component is a single edge (2 vertices); a MIS of a single edge
      // picks exactly one of its two endpoints, hence 2 options -- no need
      // to invoke the general path machinery for such a trivial shape.
      if (component.size() != 2U) {
        throw std::logic_error("2-plex complement is not a matching");
      }
      component_count = 2U;
    } else {
      // t==3: classify the component as a cycle (every vertex has
      // complement-degree exactly 2) or a path (some vertex has degree <2,
      // i.e. an endpoint) and dispatch to the matching closed-form count.
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
    // Components are independent, so their option counts multiply (see
    // enumerate_product's comment for the same cartesian-product argument).
    total = checked_multiply(total, component_count);
  }
  return total;
}

// Materializes every maximal clique of G[candidates] and invokes `callback`
// once per clique (as the full vertex-id list). This is the function the
// production Enumerator calls from try_early_termination -- see that
// function and the header's doc comment for how this fits into the overall
// BK search.
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
    // A 1-plex candidate set is already a clique in G, so it has exactly
    // one maximal clique: itself.
    callback(candidates);
    return;
  }

  const auto complement = build_complement(graph, candidates);
  std::vector<unsigned char> visited(candidates.size(), 0U);
  // Vertices forced into every output clique (complement-isolated, i.e. no
  // conflicting candidate) go straight into `fixed`; everything else is
  // grouped into per-component option lists and combined via
  // enumerate_product's cartesian product below.
  std::vector<int> fixed;
  std::vector<CliqueOptions> component_options;

  // Same DFS-over-complement-components structure as
  // count_t_plex_maximal_cliques, but building actual vertex lists instead
  // of just sizes.
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
      // A single complement edge: the two options are "pick the first
      // endpoint" or "pick the second" (no ordering/path machinery needed
      // for a 2-vertex component).
      if (component.size() != 2U) {
        throw std::logic_error("2-plex complement is not a matching");
      }
      component_options.push_back(
          {{candidates[static_cast<std::size_t>(component[0])]},
           {candidates[static_cast<std::size_t>(component[1])]}});
      continue;
    }

    // t==3: classify as cycle vs. path (see count_t_plex_maximal_cliques),
    // linearize the component into actual walk order via
    // ordered_component, then enumerate that path's/cycle's MIS options.
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

  // Combine the forced-fixed vertices with one choice per component,
  // invoking `callback` once for every combination -- i.e. once per
  // maximal clique of G[candidates].
  enumerate_product(component_options, 0U, fixed, callback);
}

} // namespace hbbmc_faithful
