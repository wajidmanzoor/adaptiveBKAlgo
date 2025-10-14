#pragma once

#include "common.h"
#include "graph.h"

// Legacy bit-based implementation (limited to 64 vertices)
void BronKerbosch(ull R, ull P, ull X);
void BKstandard(Graph& g);
void initializeAdjacencyMatrix(Graph& g);

// New dynamic implementation for any graph size
class DynamicBK {
private:
    ui n;  // number of vertices
    vector<vector<bool>> adjacency;  // adjacency matrix
    ui clique_count;
    
    void bronKerboschRecursive(vector<bool>& R, vector<bool>& P, vector<bool>& X);
    vector<bool> intersect(const vector<bool>& set1, const vector<bool>& set2);
    vector<bool> setUnion(const vector<bool>& set1, const vector<bool>& set2);
    vector<bool> setDifference(const vector<bool>& set1, const vector<bool>& set2);
    bool isEmpty(const vector<bool>& set);
    void printSet(const vector<bool>& set, const string& name);
    
public:
    DynamicBK(Graph& g);
    void findAllMaximalCliques();
    ui getCliqueCount() const { return clique_count; }
};

// Alternative implementation using adjacency lists (more memory efficient for sparse graphs)
class SparseDynamicBK {
private:
    ui n;
    vector<vector<ui>> adjacency;  // adjacency lists
    ui clique_count;
    
    void bronKerboschRecursive(vector<ui>& R, vector<ui>& P, vector<ui>& X);
    vector<ui> intersect(const vector<ui>& set1, const vector<ui>& neighbors);
    bool isEmpty(const vector<ui>& set);
    
public:
    SparseDynamicBK(Graph& g);
    void findAllMaximalCliques();
    ui getCliqueCount() const { return clique_count; }
};