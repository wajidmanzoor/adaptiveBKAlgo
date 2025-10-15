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

// Adaptive Bron-Kerbosch with Skip-Mask Tree Reordering
class AdaptiveSkipPivotBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  ui maxCliqueSize;
  
  // Skip-mask and reordering components
  vector<bool> skip_mask;     // marks vertices already covered by maximal cliques
  vector<ui> global_order;    // dynamic vertex processing order
  
  // Core algorithm methods
  void adaptiveBronKerbosch(vector<ui> &R, vector<ui> &P, vector<ui> &X, ui start_idx);
  
  // Utility methods from PivotBK
  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void applyDegreePruning(vector<ui> &P, const vector<ui> &R);
  
  // Skip-mask specific methods
  bool isClique(const vector<ui> &R) const;
  void reorderAfterClique(const vector<ui> &clique);
  vector<ui> getOrderedCandidates(const vector<ui> &P, ui start_idx) const;
  
public:
  AdaptiveSkipPivotBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};

// Simple Adaptive Enumeration with Skip-Mask Tree Reordering (No PXR/Pivot/Pruning)
class SimpleAdaptiveBK {
private:
  ui n;
  vector<vector<ui>> adjList; // adjacency lists
  ui cliqueCount;
  ui maxCliqueSize;
  
  // Skip-mask and reordering components
  vector<bool> skip_mask;     // marks vertices already covered by maximal cliques
  vector<ui> global_order;    // dynamic vertex processing order
  
  // Core enumeration method
  void enumerate(vector<ui> &R, ui start_idx);
  
  // Utility methods
  bool isConnected(ui u, ui v) const;
  bool isClique(const vector<ui> &R) const;
  void handleClique(const vector<ui> &R);
  
public:
  SimpleAdaptiveBK(Graph &g);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};
