#pragma once

#include "fast_adj_hash.h"
#include "checked_count.h"
#include "fast_clique_sink.h"

struct FastListBKTestAccess;

// Serial list kernel with adaptive ReorderSib events. Most nodes use Tomita
// pivot expansion; selected nodes use a discovered maximal clique to reorder
// the live siblings and branch on P \\ C (the single-clique sibling effect).
class FastListBK {
  friend struct FastListBKTestAccess;

private:
  static constexpr ui TINY_P_LIMIT = 4;

  struct Level {
    std::vector<ui> p;
    std::vector<ui> x;
    std::vector<ui> branch;
    std::vector<ui> processedRoots;
    std::vector<ui> witness;
  };

  const Graph &graph;
  FastAdjacencyHash adjacency;
  std::vector<ui> rank;
  std::vector<int> label;
  std::vector<Level> levels;
  ui degeneracy;
  ui minCliqueSize;
  ull cliqueCount;
  ui maxCliqueSize;
  ull checksCount;
  bool hybridReorderSibling;
  bool enableAdvancedRules;
  bool enableTailKernels;
  bool enableLocalBitset;
  ui siblingEventBudget;
  ui siblingEvents;
  ull siblingBranchesBefore;
  ull siblingBranchesAfter;
  std::vector<ui> cliqueStack;
  ull tinyKernelCalls;
  ull localBitsetHandoffs;
  ull localBitsetChecks;
  ull plex3Terminals;
  ull plex3Cliques;
  ull xDominanceRemoved;
  ull universalPForces;
  ull degreeZeroTerminals;
  ull degreeOneTerminals;
  FastCliqueSink cliqueSink;

#ifdef FASTLIST_OPPORTUNITY_PROFILE
  static constexpr size_t PROFILE_BUCKET_COUNT = 13;
  static constexpr size_t PROFILE_THRESHOLD_COUNT = 7;
  inline static constexpr std::array<ui, PROFILE_THRESHOLD_COUNT>
      PROFILE_THRESHOLDS{{3, 8, 16, 32, 64, 128, 256}};

  struct OpportunityProfile {
    std::array<ull, PROFILE_BUCKET_COUNT> stateNodes{};
    std::array<ull, PROFILE_BUCKET_COUNT> pivotItems{};
    std::array<ull, PROFILE_BUCKET_COUNT> intersectItems{};
    std::array<ull, 17 * 17> tinyPX{};

    std::array<ull, PROFILE_THRESHOLD_COUNT> fullEligible{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullEmptyCrossings{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullFrontiers{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullFrontierP{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullFrontierX{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullCells{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> fullWords{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialEligible{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialEmptyCrossings{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialFrontiers{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialFrontierP{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialFrontierX{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialCells{};
    std::array<ull, PROFILE_THRESHOLD_COUNT> partialWords{};

    ull leafMaximal = 0;
    ull leafBlockedX = 0;
    ull xUniversalPrune = 0;
    ull pCliqueSolved = 0;
    ull matchingSolved = 0;
    ull matchingOverflowFallback = 0;
    ull ordinaryBranchNodes = 0;
    ull ordinaryBranchVertices = 0;
    ull ordinaryTiny3Nodes = 0;
    ull tinyKernelNodes = 0;
    ull tinyKernelCliques = 0;

    ull pDeg0Vertices = 0;
    ull pDeg0Nodes = 0;
    ull pDeg0OutsideTiny3 = 0;
    ull pDeg1Vertices = 0;
    ull pDeg1Nodes = 0;
    ull pDeg1OutsideTiny3 = 0;
    ull universalPVertices = 0;
    ull universalPNodes = 0;
    ull universalPOutsideTiny3 = 0;

    ull plex3Nodes = 0;
    ull plex3PVertices = 0;
    ull plex3ComplementEdges = 0;
    ull plex3TopmostNodes = 0;
    ull plex3DescendantChecks = 0;
    ull plex3DescendantCliques = 0;
    ull plex3DescendantPivotItems = 0;
    ull plex3SubtreeIntersectItems = 0;

    ull pivotCandidatesX = 0;
    ull pivotCandidatesP = 0;
    ull pivotHashItemsX = 0;
    ull pivotHashItemsP = 0;
    ull pivotCsrItemsX = 0;
    ull pivotCsrItemsP = 0;
    ull pivotAbortableX = 0;
    ull pivotAbortableP = 0;
    ull pivotAbortSavedX = 0;
    ull pivotAbortSavedP = 0;
    ull naudeWidth2Nodes = 0;
    ull naudeTailCandidates = 0;
    ull naudeTailItems = 0;
    ull pImprovedPivotNodes = 0;
    ull pivotWinnerX = 0;
    ull pivotWinnerP = 0;

    ull intersectHashCalls = 0;
    ull intersectHashItems = 0;
    ull intersectCsrCalls = 0;
    ull intersectCsrItems = 0;

    ull xDomEligibleSmallNodes = 0;
    ull xDomSampledNodes = 0;
    ull xDomWarmupNodes = 0;
    ull xDomNodesWithRemoval = 0;
    ull xDomBefore = 0;
    ull xDomRemoved = 0;
    ull xDomSubsetTests = 0;

    ull graphDegree0 = 0;
    ull graphDegree1 = 0;
    ull graphDegree2 = 0;
    ull graphUndirectedEdges = 0;
    ull graphNonTriangleEdges = 0;
    ull graphNonTriangleProbes = 0;
  } profile;

  std::vector<ui> profilePAtDepth;
  std::vector<ui> profileStateAtDepth;
  bool profilePlex3Active = false;

  static size_t profileBucket(ull value);
  ull profilePivotItemsTotal() const;
  ull profileIntersectItemsTotal() const;
  void profileGraphRules();
  void profileXDominance(ui depth, const Level &level);
  void printOpportunityProfile() const;
#endif

  void buildDegeneracyOrder();
  ui neighborsInP(ui u, ui depth, const std::vector<ui> &p,
                  bool haveIncumbent, ui incumbent, bool candidateFromX);
  ui neighborsInPBaseline(ui u, ui depth,
                          const std::vector<ui> &p) const;
  bool solveTinyP(ui cliqueSize, Level &level);
  void reduceDominatedX(ui depth, Level &level);
  bool solveLowDegreeChild(ui cliqueSize, Level &child, bool &found);
  void intersectInto(ui u, ui depth, const Level &parent, Level &child);
  bool enumerateBaseline(ui depth, ui cliqueSize);
  bool enumerate(ui depth, ui cliqueSize);
  bool trySiblingEffect(ui depth, size_t nextBranch,
                        const std::vector<ui> &witness);
  bool needsCliqueStack() const {
    return hybridReorderSibling || static_cast<bool>(cliqueSink);
  }
  void emitClique(const std::vector<ui> &extension = {}) const;
  void emitComplementMatching(const std::vector<ui> &p) const;

public:
  explicit FastListBK(const Graph &g, bool hybridReorderSibling = false,
                      ui minCliqueSize = 3);
  // Installs an opt-in validation/output hook. The default empty sink keeps
  // production enumeration count-only and avoids clique materialization.
  void setCliqueSink(FastCliqueSink sink) { cliqueSink = std::move(sink); }
  void findAllMaximalCliques(const std::string &outputLabel = "FastListBK");
  ull getCliqueCount() const { return cliqueCount; }
  ui getMaxCliqueSize() const { return maxCliqueSize; }
  ull getTinyKernelCalls() const { return tinyKernelCalls; }
  ull getLocalBitsetHandoffs() const { return localBitsetHandoffs; }
  ull getPlex3Terminals() const { return plex3Terminals; }
  ull getDegreeZeroTerminals() const { return degreeZeroTerminals; }
  ull getDegreeOneTerminals() const { return degreeOneTerminals; }
};
