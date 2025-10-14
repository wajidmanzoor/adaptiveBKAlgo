#pragma once

#include "common.h"
#include "graph.h"

void BronKerbosch(ull R, ull P, ull X);
void BKstandard(Graph &g);
void initializeAdjacencyMatrix(Graph &g);

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

// Optimized Bron-Kerbosch with Pivoting and Pruning
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

// Bron-Kerbosch with Tree Reordering Strategy
class ReorderBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  set<set<ui>> foundCliques; // Store found maximal cliques

  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);
  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  
  // Subset checking functions
  bool isSubset(const set<ui> &subset, const set<ui> &superset);
  bool isSubsetOfFoundClique(const vector<ui> &candidate);
  set<ui> vectorToSet(const vector<ui> &vec);

public:
  ReorderBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
};