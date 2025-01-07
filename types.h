#pragma once

#include <cstddef>
#include <vector>

typedef enum KahyparObjectiveFunction {
  HUND_KM1,
  HUND_CUT,
  HUND_SOED,
} kahypar_objective_function_t;

typedef enum MatrixFileFormat {
  MATLAB,
  RUTHERFORD_BOEING,
  MATRIX_MARKET
} matrix_file_format_t;

struct RowColPermutation {
  std::vector<unsigned long> row_permutation;
  std::vector<unsigned long> column_permutation;
};

struct Hypergraph {
  size_t vertex_count;
  std::vector<unsigned long> hyperedges;
  std::vector<size_t> hyperedge_indices;

  size_t get_vertex_count() {
    return vertex_count;
  }

  size_t get_edge_count() {
    return hyperedge_indices.size() - 1;
  }
};