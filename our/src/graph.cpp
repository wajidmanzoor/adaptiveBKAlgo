#include "../inc/graph.h"

#include <limits>
#include <stdexcept>

namespace {

class FastIntScanner {
private:
  static constexpr size_t BUFFER_SIZE = 1 << 20;
  std::ifstream input;
  std::vector<char> buffer;
  size_t position = 0;
  size_t limit = 0;

  bool refill() {
    if (position < limit)
      return true;
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    limit = static_cast<size_t>(input.gcount());
    position = 0;
    return limit != 0;
  }

  int peek() {
    if (!refill())
      return EOF;
    return buffer[position];
  }

  int get() {
    if (!refill())
      return EOF;
    return buffer[position++];
  }

  static bool isSpace(int character) {
    return character == ' ' || character == '\n' || character == '\r' ||
           character == '\t' || character == '\v' || character == '\f';
  }

  static bool isHorizontalSpace(int character) {
    return character == ' ' || character == '\t' || character == '\v' ||
           character == '\f';
  }

  bool readDigits(ui &value) {
    unsigned long long parsed = 0;
    int character = peek();
    if (character < '0' || character > '9')
      return false;
    while (character >= '0' && character <= '9') {
      parsed = parsed * 10 + static_cast<unsigned>(character - '0');
      if (parsed > std::numeric_limits<ui>::max())
        throw std::runtime_error("graph integer exceeds uint32 range");
      get();
      character = peek();
    }
    value = static_cast<ui>(parsed);
    return true;
  }

public:
  explicit FastIntScanner(const std::string &path)
      : input(path, std::ios::in | std::ios::binary), buffer(BUFFER_SIZE) {}

  bool isOpen() const { return input.is_open(); }

  bool readUint(ui &value) {
    int character = peek();
    while (character != EOF && isSpace(character)) {
      get();
      character = peek();
    }
    return character != EOF && readDigits(value);
  }

  bool readUintOnLine(ui &value) {
    int character = peek();
    while (character != EOF && isHorizontalSpace(character)) {
      get();
      character = peek();
    }
    if (character == '\n') {
      get();
      return false;
    }
    if (character == '\r') {
      get();
      if (peek() == '\n')
        get();
      return false;
    }
    return character != EOF && readDigits(value);
  }
};

} // namespace

Graph::Graph(const std::string &path) {
  FastIntScanner scanner(path);
  if (!scanner.isOpen())
    throw std::runtime_error("cannot open graph file: " + path);
  if (!scanner.readUint(n) || !scanner.readUint(m))
    throw std::runtime_error("invalid graph header: " + path);

  if (static_cast<size_t>(m) > std::numeric_limits<size_t>::max() / 2 ||
      m > std::numeric_limits<ui>::max() / 2)
    throw std::runtime_error("graph adjacency size overflows size_t");
  const size_t expectedEntries = static_cast<size_t>(m) * 2;

  offset.assign(static_cast<size_t>(n) + 1, 0);
  degree.assign(n, 0);
  neighbors.clear();
  neighbors.reserve(expectedEntries);

  for (ui row = 0; row < n; ++row) {
    ui vertex = 0;
    if (!scanner.readUint(vertex) || vertex != row)
      throw std::runtime_error("graph rows must appear once in vertex order");

    ui neighbor = 0;
    while (scanner.readUintOnLine(neighbor)) {
      if (neighbor >= n)
        throw std::runtime_error("graph contains an out-of-range neighbor");
      if (neighbor == vertex)
        throw std::runtime_error("graph contains a self-loop");
      if (neighbors.size() == expectedEntries)
        throw std::runtime_error("graph has more adjacency entries than header");
      neighbors.push_back(neighbor);
    }

    degree[row] = static_cast<ui>(neighbors.size() - offset[row]);
    offset[row + 1] = static_cast<ui>(neighbors.size());
  }

  if (neighbors.size() != expectedEntries)
    throw std::runtime_error("graph edge count does not match header");
  ui extra = 0;
  if (scanner.readUint(extra))
    throw std::runtime_error("graph contains rows beyond declared vertex count");
}
