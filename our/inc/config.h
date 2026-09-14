#pragma once

namespace pure_config {

inline constexpr unsigned long long kDefaultBudget = 1000;
inline constexpr unsigned kHitsetCapacity = 128;
inline constexpr unsigned kCoverCollectionCutoff = 128;
// Keep the legacy PXR early-termination terminals disabled while evaluating
// the independent effect of Core Clique Removal.
inline constexpr bool kEt1Enabled = false;
inline constexpr bool kEt2Enabled = false;
inline constexpr bool kEt3Enabled = false;
inline constexpr unsigned kAdjHashThreshold = 64;
inline constexpr unsigned kSmallQCcrThreshold = 32;
inline constexpr unsigned kAdaptiveDirectQThreshold = 64;
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
