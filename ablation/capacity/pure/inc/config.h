#pragma once

namespace pure_config {

inline constexpr unsigned long long kDefaultBudget = 1000;
#if defined(PURE_HITSET_DYNAMIC)
inline constexpr bool kHitsetDynamic = true;
inline constexpr unsigned kHitsetCapacity = 0;
#else
#if !defined(PURE_HITSET_CAPACITY)
#define PURE_HITSET_CAPACITY 128
#endif
inline constexpr bool kHitsetDynamic = false;
inline constexpr unsigned kHitsetCapacity = PURE_HITSET_CAPACITY;
static_assert(kHitsetCapacity == 64 || kHitsetCapacity == 128 ||
                  kHitsetCapacity == 512 || kHitsetCapacity == 1024,
              "unsupported fixed seed-mask capacity");
#endif
// Capacity, rather than an upstream cover cap, is the varying control.
inline constexpr unsigned kCoverCollectionCutoff = 0xffffffffu;
// Keep the legacy PXR early-termination terminals disabled while evaluating
// the independent effect of Core Clique Removal.
inline constexpr bool kEt1Enabled = false;
inline constexpr bool kEt2Enabled = false;
inline constexpr bool kEt3Enabled = false;
inline constexpr unsigned kAdjHashThreshold = 256;
inline constexpr unsigned kSmallQCcrThreshold = 32;
inline constexpr unsigned kAdaptiveDirectQThreshold = 256;
inline constexpr unsigned kAdaptiveDirectWarmupRoots = 32;
inline constexpr unsigned kAdaptiveDirectMinCliquesPerRootDenominator = 2;
#if !defined(PURE_PRUNING_NORMALIZATION)
#define PURE_PRUNING_NORMALIZATION 1
#endif
#if !defined(PURE_PRUNING_SUBSUMPTION)
#define PURE_PRUNING_SUBSUMPTION 0
#endif
#if !defined(PURE_PRUNING_UNIT)
#define PURE_PRUNING_UNIT 1
#endif
#if !defined(PURE_PRUNING_USEFULNESS)
#define PURE_PRUNING_USEFULNESS 1
#endif
#if !defined(PURE_PRUNING_ANTICHAIN)
#define PURE_PRUNING_ANTICHAIN 1
#endif
#if !defined(PURE_PRUNING_FAIL_FIRST)
#define PURE_PRUNING_FAIL_FIRST 1
#endif
#if !defined(PURE_PRUNING_ZERO_COVERAGE)
#define PURE_PRUNING_ZERO_COVERAGE 1
#endif

inline constexpr bool kPruneNormalization = PURE_PRUNING_NORMALIZATION != 0;
inline constexpr bool kPruneSubsumption = PURE_PRUNING_SUBSUMPTION != 0;
inline constexpr bool kPruneUnit = PURE_PRUNING_UNIT != 0;
inline constexpr bool kPruneUsefulness = PURE_PRUNING_USEFULNESS != 0;
inline constexpr bool kPruneAntichain = PURE_PRUNING_ANTICHAIN != 0;
inline constexpr bool kPruneFailFirst = PURE_PRUNING_FAIL_FIRST != 0;
inline constexpr bool kPruneZeroCoverage = PURE_PRUNING_ZERO_COVERAGE != 0;

} // namespace pure_config
