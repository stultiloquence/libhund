#pragma once

#include <cassert>
#include <variant>
#include <vector>

#include <mtkahypartypes.h>

typedef enum class KahyparObjectiveFunction {
  KM1,
  CUT,
  SOED,
} kahypar_objective_function_t;

struct RowColPermutation {
  std::vector<unsigned long> row_permutation;
  std::vector<unsigned long> column_permutation;
};

struct BisectionConfigMtKahypar {
  double max_imbalance;
  kahypar_objective_function_t objective_function;
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

mt_kahypar_objective_t to_mt_kahypar_objective_function(
    kahypar_objective_function_t objective_function
  ) {
    switch (objective_function) {
    case KahyparObjectiveFunction::KM1: return KM1;
    case KahyparObjectiveFunction::CUT: return CUT;
    case KahyparObjectiveFunction::SOED: return SOED;
    default:
      assert(false); // Should not happen.
      return KM1;
    }
  }

  struct MultithreadingConfig {
    int number_of_threads_per_rank;
  };