#pragma once

#include "common.h"
#include "graph.h"

enum class DegOrder { ORIGINAL, ASCENDING, DESCENDING };

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
  ui checksCount;   // total vertex-set checks in bronKerboschRecursive
  ui redundancy;

  vector<vector<ui>> redendantChecks; // Track redundant vertices for pruning

  vector<vector<ui>>
      foundCliques; // Store found maximal cliques for output/debugging

  vector<ui> intersect(const vector<ui> &set1, const vector<ui> &neighbors);
  bool isEmpty(const vector<ui> &set);
  bool isConnected(ui u, ui v);
  ui choosePivot(const vector<ui> &P, const vector<ui> &X);
  void bronKerboschRecursive(vector<ui> &R, vector<ui> &P, vector<ui> &X);
  bool isPSubsetOfFoundClique(const vector<ui> &P);

public:
  PivotBK(Graph &g, DegOrder order = DegOrder::ASCENDING);

  void findAllMaximalCliques();
};

class Reorder {
private:
  ui n;
  vector<vector<ui>> adjList;
  vector<vector<ui>> adjList2;
  ui cliqueCount;
  size_t maxCliqueSize;
  ui checksCount; // total vertex-set checks in enumerate
  vector<vector<ui>> allCliques;
  vector<ui> foundLevel;
  vector<vector<ui>> cliquesByVertex;

  vector<ui> intersect(vector<ui> A, vector<ui> B);
  vector<ui> setDiff(vector<ui> A, vector<ui> B);
  vector<ui> unionSet(vector<ui> A, vector<ui> B);
  vector<ui> compliment(const vector<ui> &vector1);

  void rCall(vector<vector<ui>> mustin, vector<vector<ui>> expandTo, ui level);
  void enumerate(vector<ui> &R, vector<ui> &Q, vector<vector<ui>> &mustin,
                 vector<vector<ui>> &expandTo, ui treeIndex, ui level,
                 bool &done);

public:
  Reorder(Graph &g, DegOrder order = DegOrder::ORIGINAL);
  void findAllMaximalCliques();
  ui getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
};
