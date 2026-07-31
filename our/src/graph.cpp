#include "../inc/graph.h"

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
    return c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\v' ||
           c == '\f';
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

Graph::Graph() : n(0), m(0) {}

Graph::Graph(ui vertexCount, const std::vector<std::pair<ui, ui>> &edges)
    : n(vertexCount), m(0) {
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

Graph::Graph(std::string path) : n(0), m(0) {
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
    for (ui row = 0; row < n; row++) {
      if (!scanner.readUint(vertex))
        break;
      while (scanner.readUintOnLine(neigh)) {
        if (vertex == neigh)
          continue;
        neighbors[offset[vertex] + offset[vertex + 1]] = neigh;
        offset[vertex + 1]++;
      }
      degree[vertex] = offset[vertex + 1];
      offset[vertex + 1] += offset[vertex];
    }
  }
}
