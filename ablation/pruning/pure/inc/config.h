#pragma once

namespace pure_config {

inline constexpr unsigned long long kDefaultBudget = 1000;
inline constexpr unsigned kHitsetCapacity = 128;

#if defined(PURE_LEAN_BENCHMARK)
inline constexpr bool kDiagnosticsEnabled = false;
#else
inline constexpr bool kDiagnosticsEnabled = true;
#endif

#if defined(PURE_DISABLE_ET1)
inline constexpr bool kEt1Enabled = false;
#else
inline constexpr bool kEt1Enabled = true;
#endif
#if defined(PURE_DISABLE_ET2)
inline constexpr bool kEt2Enabled = false;
#else
inline constexpr bool kEt2Enabled = true;
#endif
#if defined(PURE_DISABLE_ET3)
inline constexpr bool kEt3Enabled = false;
#else
inline constexpr bool kEt3Enabled = true;
#endif

#if defined(PURE_ADJ_HASH_THRESHOLD)
inline constexpr unsigned kAdjHashThreshold = PURE_ADJ_HASH_THRESHOLD;
#else
inline constexpr unsigned kAdjHashThreshold = 64;
#endif

#if defined(PURE_DISABLE_PRUNING_NORMALIZATION)
inline constexpr bool kPruneNormalization = false;
#else
inline constexpr bool kPruneNormalization = true;
#endif
#if defined(PURE_DISABLE_PRUNING_SUBSUMPTION)
inline constexpr bool kPruneSubsumption = false;
#else
inline constexpr bool kPruneSubsumption = true;
#endif
#if defined(PURE_DISABLE_PRUNING_UNIT)
inline constexpr bool kPruneUnit = false;
#else
inline constexpr bool kPruneUnit = true;
#endif
#if defined(PURE_DISABLE_PRUNING_USEFULNESS)
inline constexpr bool kPruneUsefulness = false;
#else
inline constexpr bool kPruneUsefulness = true;
#endif
#if defined(PURE_DISABLE_PRUNING_ANTICHAIN)
inline constexpr bool kPruneAntichain = false;
#else
inline constexpr bool kPruneAntichain = true;
#endif
#if defined(PURE_DISABLE_PRUNING_FAIL_FIRST)
inline constexpr bool kPruneFailFirst = false;
#else
inline constexpr bool kPruneFailFirst = true;
#endif
#if defined(PURE_DISABLE_PRUNING_ZERO_COVERAGE)
inline constexpr bool kPruneZeroCoverage = false;
#else
inline constexpr bool kPruneZeroCoverage = true;
#endif

} // namespace pure_config
