#pragma once

namespace pure_config {

inline constexpr unsigned long long kDefaultBudget = 1000;
inline constexpr unsigned kHitsetCapacity = 128;
inline constexpr bool kEt1Enabled = true;
inline constexpr bool kEt2Enabled = true;
inline constexpr bool kEt3Enabled = true;
inline constexpr unsigned kAdjHashThreshold = 64;
inline constexpr unsigned kSmallQFullPxrThreshold = 4;
inline constexpr bool kPruneNormalization = true;
inline constexpr bool kPruneSubsumption = false;
inline constexpr bool kPruneUnit = true;
inline constexpr bool kPruneUsefulness = true;
inline constexpr bool kPruneAntichain = true;
inline constexpr bool kPruneFailFirst = true;
inline constexpr bool kPruneZeroCoverage = true;

} // namespace pure_config
