#pragma once

#include <cassert>
#include <variant>
#include <vector>

#include <mtkahypartypes.h>

struct RowColPermutation {
  std::vector<unsigned long> row_permutation;
  std::vector<unsigned long> column_permutation;
};

struct BisectionConfigMtKahypar {
  double max_imbalance;
};

enum BisectionConfigVariant {
  MT_KAHYPAR
};

typedef std::variant<BisectionConfigMtKahypar> BisectionConfig;

struct BreakConditionConfigRecursionDepth {
  int depth;
};

struct BreakConditionConfigBlockSize {
  int max_block_size_inclusive;
};

enum BreakConditionVariant {
  RECURSION_DEPTH,
  BLOCK_SIZE,
};

typedef std::variant<BreakConditionConfigRecursionDepth, BreakConditionConfigBlockSize> BreakConditionConfig;

  struct MultithreadingConfig {
    int number_of_threads_per_rank;
  };