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
                  kHitsetCapacity == 256 || kHitsetCapacity == 526,
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
inline constexpr bool kPruneNormalization = true;
inline constexpr bool kPruneSubsumption = false;
inline constexpr bool kPruneUnit = true;
inline constexpr bool kPruneUsefulness = true;
inline constexpr bool kPruneAntichain = true;
inline constexpr bool kPruneFailFirst = true;
inline constexpr bool kPruneZeroCoverage = true;

} // namespace pure_config
