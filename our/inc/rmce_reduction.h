#pragma once

#include "graph.h"

struct RmceReductionCounters {
  ull degree0Vertices = 0;
  ull degree1Vertices = 0;
  ull degree2Vertices = 0;
  ull nontriangleEdges = 0;
  ull directlyEmitted = 0;
  ull discardedUnderThreshold = 0;
  ull duplicateDirectOutputs = 0;
};

struct RmceReductionResult {
  Graph graph;
  vector<ui> residualToOriginal;
  ull directlyEmittedCount = 0;
  size_t maximumCliqueSize = 0;
  vector<vector<ui>> directlyEmittedCliques;
  RmceReductionCounters counters;
};

RmceReductionResult applyRmceReduction(const Graph &graph, ui minCliqueSize,
                                       bool collectCliqueIdentities);
