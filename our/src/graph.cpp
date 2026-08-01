#include "../inc/graph.h"
#include <numeric>

namespace {
class FastIntScanner {
private:
  static constexpr size_t BUFFER_SIZE = 1 << 20;
  std::ifstream in;
  std::vector<char> buffer;
  size_t pos;
  size_t limit;

  bool refill() {
    if (pos < limit)
      return true;
    in.read(buffer.data(), buffer.size());
    limit = static_cast<size_t>(in.gcount());
    pos = 0;
    return limit != 0;
  }

  int peek() {
    if (!refill())
      return EOF;
    return buffer[pos];
  }

  int get() {
    if (!refill())
      return EOF;
    return buffer[pos++];
  }

  static bool isSpace(int c) {
    return c == ' ' || c == '\n' || c == '\r' || c == '\t' ||
           c == '\v' || c == '\f';
  }

  static bool isHorizontalSpace(int c) {
    return c == ' ' || c == '\t' || c == '\v' || c == '\f';
  }

public:
  explicit FastIntScanner(const std::string &path)
      : in(path, std::ios::in | std::ios::binary), buffer(BUFFER_SIZE), pos(0),
        limit(0) {}

  bool isOpen() const { return in.is_open(); }

  bool readUint(ui &value) {
    int c = peek();
    while (c != EOF && isSpace(c)) {
      get();
      c = peek();
    }
    if (c == EOF)
      return false;
    if (c < '0' || c > '9')
      return false;

    value = 0;
    while (c >= '0' && c <= '9') {
      value = value * 10 + static_cast<ui>(c - '0');
      get();
      c = peek();
    }
    return true;
  }

  bool readUintOnLine(ui &value) {
    int c = peek();
    while (c != EOF && isHorizontalSpace(c)) {
      get();
      c = peek();
    }
    if (c == EOF)
      return false;
    if (c == '\n') {
      get();
      return false;
    }
    if (c == '\r') {
      get();
      if (peek() == '\n')
        get();
      return false;
    }
    if (c < '0' || c > '9') {
      get();
      return false;
    }

    value = 0;
    while (c >= '0' && c <= '9') {
      value = value * 10 + static_cast<ui>(c - '0');
      get();
      c = peek();
    }
    return true;
  }
};
} // namespace

Graph::Graph() : n(0), m(0), kmax(0), adjacencySorted(false) {}

Graph::Graph(ui vertexCount,
             const std::vector<std::pair<ui, ui>> &edges)
    : n(vertexCount), m(0), kmax(0), adjacencySorted(true) {
  std::vector<std::vector<ui>> rows(n);
  for (const auto &[u, v] : edges) {
    if (u >= n || v >= n || u == v)
      continue;
    rows[u].push_back(v);
    rows[v].push_back(u);
  }

  offset.assign(n + 1, 0);
  degree.assign(n, 0);
  for (ui v = 0; v < n; v++) {
    auto &row = rows[v];
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    degree[v] = static_cast<ui>(row.size());
    offset[v + 1] = offset[v] + degree[v];
  }

  neighbors.resize(offset[n]);
  for (ui v = 0; v < n; v++)
    std::copy(rows[v].begin(), rows[v].end(), neighbors.begin() + offset[v]);
  m = static_cast<ui>(neighbors.size() / 2);
}

Graph::Graph(std::string path) : n(0), m(0), kmax(0), adjacencySorted(false) {
  FastIntScanner scanner(path);

  if (!scanner.isOpen()) {
    std::cout << "Graph file Open Failed " << std::endl;
    exit(1);
  } else {
    scanner.readUint(n);
    scanner.readUint(m);

    offset.resize(n + 1, 0);
    neighbors.resize(2 * m);
    degree.resize(n, 0);
    ui vertex, neigh;
    bool rowsSorted = true;
    for (ui row = 0; row < n; row++) {
      if (!scanner.readUint(vertex))
        break;
      bool haveLastNeigh = false;
      ui lastNeigh = 0;
      while (scanner.readUintOnLine(neigh)) {
        if (vertex == neigh)
          continue;
        if (haveLastNeigh && neigh < lastNeigh)
          rowsSorted = false;
        lastNeigh = neigh;
        haveLastNeigh = true;
        neighbors[offset[vertex] + offset[vertex + 1]] = neigh;
        offset[vertex + 1]++;
      }
      degree[vertex] = offset[vertex + 1];
      offset[vertex + 1] += offset[vertex];
    }
    adjacencySorted = rowsSorted;
  }

  if (debug) {
    std::cout << "n=" << n << ", m=" << m << std::endl;

    cout << "ofsset ";
    for (ui i = 0; i <= n; i++) {
      cout << offset[i] << " ";
    }
    cout << endl;
    cout << "degree ";
    for (ui i = 0; i < n; i++) {
      cout << degree[i] << " ";
    }
    cout << endl;
    cout << "neighbors ";
    for (ui i = 0; i < 2 * m; i++) {
      cout << neighbors[i] << " ";
    }
    cout << endl;
  }
}

void Graph::sortAdjacency() {
  if (adjacencySorted)
    return;

  for (ui v = 0; v < n; v++) {
    auto first = neighbors.begin() + offset[v];
    auto last = neighbors.begin() + offset[v + 1];
    if (!std::is_sorted(first, last))
      std::sort(first, last);
  }
  adjacencySorted = true;
}

void Graph::getListingOrder(std::vector<ui> &arr) {
  /* Rettrun an array with each index storing the listing order based on the
  core value. Listing order is a unique number. high core values get low listing
  order*/
  corePeelSequence.resize(n);
  coreDecompose(corePeelSequence);

  for (size_t i = 0; i < n; ++i) {
    arr[corePeelSequence[i]] = i + 1;
  }
}

void Graph::coreDecompose(std::vector<ui> &arr) {
  /* Peeling algorithm to find the core values of each vertex.
     Returns the peeling sequence i.e. verticies in increaseing order of core
     values. */
  core.resize(n);
  if (n == 0) {
    kmax = 0;
    return;
  }

  int maxDegree = *std::max_element(degree.begin(), degree.end());
  if (debug) std::cout << "Max Degree=" << maxDegree << std::endl;

  // Initialize bins
  std::vector<ui> bins(maxDegree + 1, 0);
  for (ui deg : degree) {
    bins[deg]++;
  }

  // Compute bin positions
  std::vector<int> bin_positions(maxDegree + 1, 0);
  std::partial_sum(bins.begin(), bins.end(), bin_positions.begin());

  // Initialize position and sortedVertex arrays
  std::vector<ui> position(n);
  std::vector<ui> sortedVertex(n);

  for (ui v = 0; v < n; v++) {
    position[v] = --bin_positions[degree[v]]; // Assign position
    sortedVertex[position[v]] = v;            // Place vertex in sorted list
  }

  // Perform core decomposition
  for (int i = 0; i < n; i++) {
    ui v = sortedVertex[i];
    core[v] = degree[v]; // Assign core value
    arr[n - i - 1] = v;  // Assign peel sequence

    // Update degrees of neighbors
    for (int j = offset[v]; j < offset[v + 1]; j++) {
      ui u = neighbors[j];
      if (degree[u] > degree[v]) {
        ui du = degree[u];
        ui pu = position[u];
        ui pw = bin_positions[du];
        ui w = sortedVertex[pw];

        if (u != w) {
          position[u] = pw;
          sortedVertex[pu] = w;
          position[w] = pu;
          sortedVertex[pw] = u;
        }

        bin_positions[du]++;
        degree[u]--;
      }
    }
  }

  kmax = core[0]; // Initialize with first element
  for (ui val : core) {
    if (val > kmax)
      kmax = val;
  }
}
