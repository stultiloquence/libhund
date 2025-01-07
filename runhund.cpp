#include <cstdio>
#include <fstream>
#include <memory>
#include <iostream>
#include <vector>
#include <thread>

#include "context.h"

#include "types.h"
#include "util.h"

struct triplet_matrix {
    int64_t nrows = 0, ncols = 0;
    std::vector<int64_t> rows, cols;
    std::vector<double> vals;       // or int64_t, float, std::complex<double>, etc.
} mat;

int main(int argc, char const *argv[]) {
	std::vector<size_t> example_hyperedge_indices { 0, 3, 5, 8, 11, 15, 17, 19, 22, 24 };
	std::vector<unsigned long> example_hyperedges { 0, 4, 7, 1, 3, 2, 4, 7, 3, 5, 6, 1, 2, 3, 4, 1, 5, 5, 6, 2, 3, 7, 0, 7 };

	auto ctx = Context(
		(BisectionConfigMtKahypar) {
			.max_imbalance = 0.03,
			.objective_function = HUND_KM1,
		},
		(BreakConditionConfigBlockSize) {
			.block_size = 5,
		},
		(HypergraphConfigFromFile) {
			.file_format = MATRIX_MARKET,
			.file_path = "1138_bus.mtx",
		}
	);

	auto result = ctx.run();

	print_vector(result.row_permutation);
	print_vector(result.column_permutation);
}
