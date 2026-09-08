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

#if defined(PURE_LEAN_BENCHMARK)
inline void addSearchStatesOrThrow(ull &, ull) noexcept {}

inline void incrementSearchStateOrThrow(ull &) noexcept {}
#else
inline void addSearchStatesOrThrow(ull &count, ull increment) {
  ull sum = 0;
  if (!tryAddUll(count, increment, sum))
    throw std::overflow_error(
        "search-state count exceeds the uint64_t range");
  count = sum;
}

inline void incrementSearchStateOrThrow(ull &count) {
  addSearchStatesOrThrow(count, 1);
}
#endif
