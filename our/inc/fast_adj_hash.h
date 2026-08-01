#pragma once

#include "graph.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

#if defined(__SSE2__) &&                                                   \
    (defined(__x86_64__) || defined(__i386__) || defined(_M_X64) ||       \
     defined(_M_IX86))
#include <emmintrin.h>
#define FAST_ADJACENCY_HASH_SSE2 1
#else
#define FAST_ADJACENCY_HASH_SSE2 0
#endif

// Static adjacency membership adapted from HBBMC's CuckooHash. Each key has
// two complementary four-slot buckets. Rows are built once, stored in one
// contiguous array, and queried with two aligned SSE2 loads when available.
// The scalar path uses the identical layout on ARM and other targets.
class FastAdjacencyHash {
private:
  static constexpr ui EMPTY = std::numeric_limits<ui>::max();
  static constexpr ui BUCKET_WIDTH = 4;

  struct alignas(16) Bucket {
    ui value[BUCKET_WIDTH];

    Bucket() : value{EMPTY, EMPTY, EMPTY, EMPTY} {}
  };

  struct Row {
    size_t offset = 0; // Bucket offset, not byte/element offset.
    ui mask = 0;
  };

  std::vector<Row> rows;
  std::vector<Bucket> buckets;

  static bool bucketContains(const Bucket &bucket, ui value) {
    for (ui i = 0; i < BUCKET_WIDTH; ++i) {
      if (bucket.value[i] == value)
        return true;
    }
    return false;
  }

  // Fold high vertex-ID bits into the bucket index. Raw low-bit indexing can
  // otherwise allocate enormous rows for a handful of regularly spaced IDs.
  static ui mixedValue(ui value) {
    return (value * 0x9e3779b1U) ^ (value >> 16);
  }

  static ui firstSlot(ui value, ui mask) {
    return mixedValue(value) & mask;
  }

  static ui secondSlot(ui value, ui mask) {
    return (~mixedValue(value)) & mask;
  }

  static void doubleBucketCount(ui &bucketCount) {
    if (bucketCount > std::numeric_limits<ui>::max() / 2)
      throw std::overflow_error("adjacency hash bucket count overflow");
    bucketCount <<= 1;
  }

  static bool contains(const std::vector<Bucket> &row, ui mask, ui value) {
    return bucketContains(row[firstSlot(value, mask)], value) ||
           bucketContains(row[secondSlot(value, mask)], value);
  }

  static bool insert(std::vector<Bucket> &row, ui mask, ui value) {
    if (contains(row, mask, value))
      return true;

    ui first = firstSlot(value, mask);
    for (ui i = 0; i < BUCKET_WIDTH; ++i) {
      if (row[first].value[i] == EMPTY) {
        row[first].value[i] = value;
        return true;
      }
    }

    ui second = secondSlot(value, mask);
    for (ui i = 0; i < BUCKET_WIDTH; ++i) {
      if (row[second].value[i] == EMPTY) {
        row[second].value[i] = value;
        return true;
      }
    }

    // Same deterministic eviction policy as compare/HBBMC/src/hash.hpp.
    bool useFirst = true;
    ui current = value;
    for (ui attempt = 0; attempt < mask; ++attempt) {
      const ui slot =
          useFirst ? firstSlot(current, mask) : secondSlot(current, mask);
      Bucket &bucket = row[slot];
      ui empty = 0;
      while (empty < BUCKET_WIDTH && bucket.value[empty] != EMPTY)
        ++empty;

      ui replaced;
      if (empty == BUCKET_WIDTH) {
        replaced = bucket.value[0];
        for (ui i = 1; i < BUCKET_WIDTH; ++i)
          bucket.value[i - 1] = bucket.value[i];
        bucket.value[BUCKET_WIDTH - 1] = current;
      } else {
        replaced = EMPTY;
        bucket.value[empty] = current;
      }
      useFirst = slot == secondSlot(replaced, mask);
      current = replaced;
      if (current == EMPTY)
        return true;
    }
    return false;
  }

public:
  explicit FastAdjacencyHash(const Graph &g) : rows(g.n) {
    static_assert(sizeof(ui) == sizeof(std::uint32_t),
                  "FastAdjacencyHash requires 32-bit vertex IDs");
    static_assert(sizeof(Bucket) == 16 && alignof(Bucket) >= 16,
                  "FastAdjacencyHash buckets must be 16-byte SIMD lanes");

    size_t estimatedBuckets = 0;
    for (ui degree : g.degree) {
      ui bucketCount = 2;
      while (static_cast<ull>(degree) + 1 >=
             (static_cast<ull>(bucketCount) - 1) * BUCKET_WIDTH)
        doubleBucketCount(bucketCount);
      estimatedBuckets += bucketCount;
    }
    buckets.reserve(estimatedBuckets);

    for (ui u = 0; u < g.n; ++u) {
      ui bucketCount = 2;
      while (static_cast<ull>(g.degree[u]) + 1 >=
             (static_cast<ull>(bucketCount) - 1) * BUCKET_WIDTH)
        doubleBucketCount(bucketCount);

      std::vector<Bucket> row;
      for (;;) {
        row.assign(bucketCount, Bucket{});
        bool built = true;
        for (ui at = g.offset[u]; at < g.offset[u + 1]; ++at) {
          if (!insert(row, bucketCount - 1, g.neighbors[at])) {
            built = false;
            break;
          }
        }
        if (built)
          break;
        doubleBucketCount(bucketCount);
      }

      rows[u].offset = buckets.size();
      rows[u].mask = bucketCount - 1;
      buckets.insert(buckets.end(), row.begin(), row.end());
    }
  }

  bool contains(ui u, ui v) const {
    const Row &row = rows[u];
    const Bucket &first = buckets[row.offset + firstSlot(v, row.mask)];
    const Bucket &second = buckets[row.offset + secondSlot(v, row.mask)];
#if FAST_ADJACENCY_HASH_SSE2
    const __m128i needle = _mm_set1_epi32(static_cast<int>(v));
    const __m128i firstValues =
        _mm_load_si128(reinterpret_cast<const __m128i *>(first.value));
    const __m128i secondValues =
        _mm_load_si128(reinterpret_cast<const __m128i *>(second.value));
    const __m128i matches =
        _mm_or_si128(_mm_cmpeq_epi32(needle, firstValues),
                     _mm_cmpeq_epi32(needle, secondValues));
    return _mm_movemask_epi8(matches) != 0;
#else
    return bucketContains(first, v) || bucketContains(second, v);
#endif
  }
};

#undef FAST_ADJACENCY_HASH_SSE2
