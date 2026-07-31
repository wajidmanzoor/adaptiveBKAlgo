#pragma once

#include "common.h"

#include <limits>
#include <stdexcept>

static_assert(std::numeric_limits<ull>::digits == 64,
              "clique counters require a 64-bit ull type");

inline bool tryAddUll(ull left, ull right, ull &result) noexcept {
  if (right > std::numeric_limits<ull>::max() - left)
    return false;
  result = left + right;
  return true;
}

inline bool tryMultiplyUll(ull left, ull right, ull &result) noexcept {
  if (left != 0 && right > std::numeric_limits<ull>::max() / left)
    return false;
  result = left * right;
  return true;
}

inline void addCliqueCountOrThrow(ull &count, ull increment) {
  ull sum = 0;
  if (!tryAddUll(count, increment, sum))
    throw std::overflow_error(
        "maximal-clique count exceeds the uint64_t output range");
  count = sum;
}

inline void incrementSearchStateOrThrow(ull &count) {
  if (count == std::numeric_limits<ull>::max())
    throw std::overflow_error(
        "recursive search-state count exceeds the uint64_t range");
  ++count;
}
