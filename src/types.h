#pragma once

#include <cassert>
#include <variant>
#include <vector>

struct RowColPermutation {
  std::vector<unsigned long> row_permutation;
  std::vector<unsigned long> column_permutation;
};

struct BreakConditionConfigRecursionDepth {
  int depth;
};

enum BisectionConfigVariant { MT_KAHYPAR };

struct BreakConditionConfigBlockSize {
  int max_block_size_inclusive;
};

enum BreakConditionVariant {
  RECURSION_DEPTH,
  BLOCK_SIZE,
};

typedef std::variant<BreakConditionConfigRecursionDepth, BreakConditionConfigBlockSize> BreakConditionConfig;