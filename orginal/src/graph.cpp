#include "../inc/graph.h"

Graph::Graph() {
  // Default constructor implementation
}

Graph::Graph(std::string path) {
  std::ifstream inputFile(path, std::ios::in);

  if (!inputFile.is_open()) {
    std::cout << "Graph file Open Failed " << std::endl;
    exit(1);
  }
  
  std::string line;
  std::getline(inputFile, line);
  std::istringstream iss(line);
  iss >> n >> m;

  offset.resize(n + 1, 0);
  neighbors.resize(2 * m);
  degree.resize(n, 0);
  
  // First pass: count degrees
  int vertex, neigh;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    iss >> vertex;
    while (iss >> neigh) {
      if (vertex != neigh) {
        degree[vertex]++;
      }
    }
  }
  
  // Compute offsets
  for (ui i = 1; i <= n; i++) {
    offset[i] = offset[i-1] + degree[i-1];
  }
  
  // Reset degrees for second pass
  std::fill(degree.begin(), degree.end(), 0);
  
  // Second pass: fill neighbors
  inputFile.clear();
  inputFile.seekg(0);
  std::getline(inputFile, line); // Skip first line
  
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    iss >> vertex;
    while (iss >> neigh) {
      if (vertex != neigh) {
        neighbors[offset[vertex] + degree[vertex]] = neigh;
        degree[vertex]++;
      }
    }
  }

  inputFile.close();
  std::cout << "n=" << n << ", m=" << m << std::endl;
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
  int maxDegree = *std::max_element(degree.begin(), degree.end());
  std::cout << "Max Degree=" << maxDegree << std::endl;

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

void Graph::convertToBitmap(ull* pik, byte& N) {
  N = n;
  
  // Initialize pik array
  for (ui i = 0; i < n; i++) {
    pik[i] = 0;
  }
  
  // Convert adjacency list to bitmap representation
  for (ui i = 0; i < n; i++) {
    ull mask = 1ULL << (n - 1 - i);  // Bit for vertex i
    for (ui j = offset[i]; j < offset[i + 1]; j++) {
      ui neighbor = neighbors[j];
      pik[i] |= (1ULL << (n - 1 - neighbor));  // Set bit for neighbor
    }
  }
}