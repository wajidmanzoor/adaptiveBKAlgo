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
  PivotBK(Graph &g);

  void findAllMaximalCliques();
};

class ReorderBK2 {
private:
  Graph graph;
  ui n;
  vector<vector<ui>> adjList;  // adjacency lists
  vector<vector<ui>> adjList2; // adjacency lists with only greater vertices
  ui cliqueCount;
  size_t maxCliqueSize;
  vector<ui> order; // ordering of vertices
  vector<ui> old2new;
  vector<ui> new2old;
  vector<vector<ui>> foundCliques;
  // vector<bool> visited; // tracks which vertices have been starting points
  //  vector<bool> status;  // tracks status of vertices
  // bool canExtend(const vector<ui> &R, ui vertex) const;
  // bool isConnected(ui u, ui v) const;

  // Set operations
  vector<ui> intersect(vector<ui> vector1, vector<ui> vector2);
  vector<ui> setDifference(vector<ui> A, vector<ui> B);
  vector<ui> unionSet(vector<ui> vector1, vector<ui> vector2);

  // Calls expansion of Trees
  void rCall(vector<vector<ui>> &mustin, vector<vector<ui>> &expandTo, ui level,
             ui enlevel);

  // Enumerates a tree untill a Maximal clique is found
  void enemurate(vector<ui> &R, vector<ui> &Q, vector<vector<ui>> &mustin,
                 vector<vector<ui>> &expandTo, bool &moveToNext, bool &flag,
                 ui index, ui level, ui enlevel);

public:
  ReorderBK2(Graph &g, bool sortMode = true);
  // Finds all Maximal Cliques
  void findAllMaximalCliques();
  // Returns number of Maximal cliques found
  ui getCliqueCount() const { return cliqueCount; }
  // Returns size of maximum maximal clique found
  ui getMaxCliqueSize() const { return maxCliqueSize; }
  vector<ui> applyPreviousCliqueEffects(const vector<ui> &tree,
                                        const vector<ui> &Q);
};
