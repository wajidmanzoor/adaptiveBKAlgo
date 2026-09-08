#pragma once

#include "common.h"

class Graph {
public:
  ui n = 0;
  ui m = 0;

  std::vector<ui> offset;
  std::vector<ui> neighbors;
  std::vector<ui> degree;

  explicit Graph(const std::string &path);
};
