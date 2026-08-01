#pragma once

#include "common.h"

class Graph {
public:
  ui n;
  ui m;
  ui kmax;

  std::vector<ui> offset;
  std::vector<ui> neighbors;
  std::vector<ui> degree;
  std::vector<ui> core;
  std::vector<ui> corePeelSequence;
  std::string filePath;
  bool adjacencySorted;

public:
  Graph();
  Graph(std::string path);
  Graph(ui vertexCount, const std::vector<std::pair<ui, ui>> &edges);
  void sortAdjacency();
  void getListingOrder(std::vector<ui> &arr);
  void coreDecompose(std::vector<ui> &arr);
};
