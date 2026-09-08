#include "../inc/graph.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace hbbmc_faithful {

namespace {

// Combines an unordered vertex pair into one 64-bit hash key: smaller id in
// the high 32 bits, larger id in the low 32 bits. Swapping first makes the
// key symmetric, i.e. key(u,v) == key(v,u), which is required since Edge
// storage always normalizes to u < v.
std::uint64_t local_edge_key(int u, int v) {
  if (u > v) {
    std::swap(u, v);
  }
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32U) |
         static_cast<std::uint32_t>(v);
}

} // namespace

Graph::Graph(int vertex_count, const std::vector<std::pair<int, int>> &edges) {
  std::vector<std::int64_t> labels(static_cast<std::size_t>(vertex_count));
  for (int v = 0; v < vertex_count; ++v) {
    labels[static_cast<std::size_t>(v)] = v;
  }
  build(vertex_count, edges, std::move(labels));
}

Graph::Graph(std::vector<std::int64_t> original_labels,
             const std::vector<std::pair<int, int>> &edges) {
  const int vertex_count = static_cast<int>(original_labels.size());
  build(vertex_count, edges, std::move(original_labels));
}

Graph Graph::read_edge_list(const std::string &path, int declared_vertices) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open graph file: " + path);
  }

  // First pass: collect every "u v" edge with its *original* (possibly huge
  // or sparse) labels, plus the flat list of labels seen, so we can later
  // decide how to map labels onto dense ids.
  std::vector<std::pair<std::int64_t, std::int64_t>> labeled_edges;
  std::vector<std::int64_t> labels;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#' ||
        line[first] == '%') {
      continue;
    }
    std::istringstream row(line);
    std::int64_t u = 0;
    std::int64_t v = 0;
    if (!(row >> u >> v)) {
      throw std::runtime_error("invalid edge at " + path + ":" +
                               std::to_string(line_number));
    }
    if (u == v) {
      continue;
    }
    labeled_edges.emplace_back(u, v);
    labels.push_back(u);
    labels.push_back(v);
  }

  Graph graph;
  if (declared_vertices >= 0) {
    // --num-vertices mode: the caller asserts the label space is exactly
    // 0..declared_vertices-1, so labels ARE dense ids already. This is the
    // only path that can produce isolated vertices (a label that never
    // appears in an edge still gets a slot), which matters because HBBMC
    // reports isolated vertices as maximal 1-cliques (see enumerator.cpp).
    for (const auto &[u, v] : labeled_edges) {
      if (u < 0 || v < 0 || u >= declared_vertices || v >= declared_vertices) {
        throw std::runtime_error(
            "--num-vertices requires every input label to be in [0,N)");
      }
    }
    std::vector<std::pair<int, int>> dense_edges;
    dense_edges.reserve(labeled_edges.size());
    for (const auto &[u, v] : labeled_edges) {
      dense_edges.emplace_back(static_cast<int>(u), static_cast<int>(v));
    }
    std::vector<std::int64_t> dense_labels(
        static_cast<std::size_t>(declared_vertices));
    for (int v = 0; v < declared_vertices; ++v) {
      dense_labels[static_cast<std::size_t>(v)] = v;
    }
    graph.build(declared_vertices, dense_edges, std::move(dense_labels));
    return graph;
  }

  // Default mode: only labels that actually occur in some edge become
  // vertices. Sort+unique the collected labels to get the vertex universe,
  // then remap every edge endpoint through a label->dense-id hash map.
  std::sort(labels.begin(), labels.end());
  labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
  if (labels.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error(
        "graph has too many vertices for 32-bit dense ids");
  }
  std::unordered_map<std::int64_t, int> dense;
  dense.reserve(labels.size() * 2U + 1U);
  for (std::size_t i = 0; i < labels.size(); ++i) {
    dense.emplace(labels[i], static_cast<int>(i));
  }
  std::vector<std::pair<int, int>> dense_edges;
  dense_edges.reserve(labeled_edges.size());
  for (const auto &[u, v] : labeled_edges) {
    dense_edges.emplace_back(dense.at(u), dense.at(v));
  }
  const int vertex_count = static_cast<int>(labels.size());
  graph.build(vertex_count, dense_edges, std::move(labels));
  return graph;
}

void Graph::build(int vertex_count,
                  const std::vector<std::pair<int, int>> &input_edges,
                  std::vector<std::int64_t> original_labels) {
  if (vertex_count < 0) {
    throw std::invalid_argument("negative vertex count");
  }
  if (original_labels.size() != static_cast<std::size_t>(vertex_count)) {
    throw std::invalid_argument("label count does not match vertex count");
  }

  adjacency_.assign(static_cast<std::size_t>(vertex_count), {});
  edges_.clear();
  edge_ids_.clear();
  original_labels_ = std::move(original_labels);

  // Normalize every edge to (min, max), then sort + dedup so self-loops and
  // repeated edges from the input collapse into one canonical edge list.
  std::vector<std::pair<int, int>> edges;
  edges.reserve(input_edges.size());
  for (auto [u, v] : input_edges) {
    if (u < 0 || v < 0 || u >= vertex_count || v >= vertex_count) {
      throw std::invalid_argument("edge endpoint outside vertex range");
    }
    if (u == v) {
      continue;
    }
    if (u > v) {
      std::swap(u, v);
    }
    edges.emplace_back(u, v);
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  // Assign dense edge ids by walking the sorted (u,v) list in order, and
  // append each endpoint to the other's adjacency list as we go.
  //
  // Why this leaves every adjacency_[x] sorted without an explicit sort
  // pass: for a fixed vertex x, edges (a,x) with a<x are only encountered
  // while the outer loop's "u" is at a, and u only ever increases across
  // the loop, so those a's are appended to adjacency_[x] in increasing
  // order. When the loop's u reaches x, all edges (x,b) with b>x are
  // processed contiguously with b already increasing (they were sorted by
  // (u,v)) and are appended after every a<x. So adjacency_[x] ends up as
  // "smaller neighbors ascending" followed by "larger neighbors ascending"
  // — i.e. fully sorted. That sortedness is what lets adjacent() below use
  // binary_search, and lets the rest of the codebase treat neighbor lists
  // as sorted sets for set_intersection/set_difference.
  edges_.reserve(edges.size());
  edge_ids_.reserve(edges.size() * 2U + 1U);
  for (const auto &[u, v] : edges) {
    const int id = static_cast<int>(edges_.size());
    edges_.push_back({u, v});
    edge_ids_.emplace(local_edge_key(u, v), id);
    adjacency_[static_cast<std::size_t>(u)].push_back(v);
    adjacency_[static_cast<std::size_t>(v)].push_back(u);
  }
}

std::uint64_t Graph::edge_key(int u, int v) { return local_edge_key(u, v); }

bool Graph::adjacent(int u, int v) const {
  if (u == v || u < 0 || v < 0 || u >= vertex_count() || v >= vertex_count()) {
    return false;
  }
  const auto &row = adjacency_[static_cast<std::size_t>(u)];
  return std::binary_search(row.begin(), row.end(), v);
}

int Graph::edge_id(int u, int v) const {
  if (u == v) {
    return -1;
  }
  const auto it = edge_ids_.find(edge_key(u, v));
  return it == edge_ids_.end() ? -1 : it->second;
}

// Relies on both adjacency lists being sorted (see build()) to intersect
// them in linear time via std::set_intersection instead of a hash lookup
// per candidate.
std::vector<int> Graph::common_neighbors(int u, int v) const {
  std::vector<int> result;
  const auto &a = adjacency_.at(static_cast<std::size_t>(u));
  const auto &b = adjacency_.at(static_cast<std::size_t>(v));
  result.reserve(std::min(a.size(), b.size()));
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::back_inserter(result));
  return result;
}

} // namespace hbbmc_faithful
