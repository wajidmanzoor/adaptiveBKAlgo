#pragma once

#include "common.h"
#include "graph.h"

// Adjacency Matrix based Bron-Kerbosch
class AdjMatBK {
private:
  ui n;
  vector<vector<bool>> adjMatrix;
  ui cliqueCount;
  void bronKerboschRecursive(vector<bool> &R, vector<bool> &P, vector<bool> &X);
  vector<bool> intersect(const vector<bool> &set1, const vector<bool> &set2);
  vector<bool> setUnion(const vector<bool> &set1, const vector<bool> &set2);
  vector<bool> setDifference(const vector<bool> &set1,
                             const vector<bool> &set2);
  bool isEmpty(const vector<bool> &set);
  void printSet(const vector<bool> &set, const string &name);

public:
  AdjMatBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
};
// Adjacency List based Bron-Kerbosch
class AdjListBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;

  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);
  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);

public:
  AdjListBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
};

// Optimized Adjacency List based Bron-Kerbosch with Pivoting and Pruning
class PivotBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  ui maxCliqueSize; // for neighborhood size pruning

  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);
  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);

  // Pivoting strategy
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);

  // Pruning rules
  bool isRedundantX(ui x, const vector<ui> &P);
  bool isConnected(ui u, ui v);
  void applyDegreePruning(vector<ui> &P, const vector<ui> &R);
  void applyNeighborhoodPruning(vector<ui> &P, const vector<ui> &R);

public:
  PivotBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
};

class ReorderBK2 {
private:
  Graph graph;
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  size_t maxCliqueSize;
  vector<bool> visited; // tracks which vertices have been starting points
  vector<bool> status;  // tracks status of vertices
  bool canExtend(const vector<ui> &R, ui vertex) const;
  bool isConnected(ui u, ui v) const;
  void rCall(vector<ui> &expandFrom, vector<ui> &expandTo);
  vector<ui> intersect(vector<ui> vector1, vector<ui> vector2);

public:
  ReorderBK2(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};

class ReorderBK {
private:
  // ----- Graph and adjacency -----
  Graph graph;
  ui n;
  vector<vector<ui>> adjList;

  // ----- Global Clique Stats -----
  ui cliqueCount;
  ui maxCliqueSize;

  // ----- Core of Recursive Re-Ordering Strategy -----

  // Global vertex ordering (critical)
  vector<ui> order; // order[i] = the vertex at position i
  vector<ui> pos;   // pos[v] = index in order[]

  // Boundary between "frontier" and "tail"
  int tail_start; // vertices in positions >= tail_start are in the tail block

  // ----- Internal Methods -----
  bool isConnected(ui u, ui v) const; // adjacency check
  bool canExtend(const vector<ui> &R, ui v) const;

  // Move vertex v to tail region once its branch is done
  void shiftToTail(ui v);

  // Exact recursive enumeration
  void search(vector<ui> &R, int startIdx);

public:
  ReorderBK(Graph &g);
  void run();

  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};
