#pragma once

#include <cassert>
#include <variant>
#include <vector>

#include <mtkahypartypes.h>

typedef enum KahyparObjectiveFunction {
  HUND_KM1,
  HUND_CUT,
  HUND_SOED,
} kahypar_objective_function_t;

struct RowColPermutation {
  std::vector<unsigned long> row_permutation;
  std::vector<unsigned long> column_permutation;
};

struct BisectionConfigMtKahypar {
  double max_imbalance;
  kahypar_objective_function_t objective_function;
};

typedef std::variant<BisectionConfigMtKahypar> BisectionConfig;

struct BreakConditionConfigRecursionDepth {
  int depth;
};

struct BreakConditionConfigBlockSize {
  int max_block_size_inclusive;
};

typedef std::variant<BreakConditionConfigRecursionDepth, BreakConditionConfigBlockSize> BreakConditionConfig;

mt_kahypar_objective_t to_mt_kahypar_objective_function(
    kahypar_objective_function_t objective_function
  ) {
    switch (objective_function) {
    case HUND_KM1: return KM1;
    case HUND_CUT: return CUT;
    case HUND_SOED: return SOED;
    default:
      assert(false); // Should not happen.
      return KM1;
    }
  }