#pragma once

#include "common.h"

class Graph {
public:
  ui n;
  ui m;

  std::vector<ui> offset;
  std::vector<ui> neighbors;
  std::vector<ui> degree;

public:
  Graph();
  explicit Graph(std::string path);
  Graph(ui vertexCount, const std::vector<std::pair<ui, ui>> &edges);
};
