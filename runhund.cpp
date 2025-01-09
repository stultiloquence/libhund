#include <cstdio>
#include <vector>

#include <mpi.h>

#include "context.h"
#include "types.h"
#include "util.h"

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	std::vector<size_t> example_hyperedge_indices { 0, 3, 5, 8, 11, 15, 17, 19, 22, 24 };
	std::vector<unsigned long> example_hyperedges { 0, 4, 7, 1, 3, 2, 4, 7, 3, 5, 6, 1, 2, 3, 4, 1, 5, 5, 6, 2, 3, 7, 0, 7 };

	auto ctx = Context(
		(BisectionConfigMtKahypar) {
			.max_imbalance = 0.03,
			.objective_function = HUND_KM1,
		},
		(BreakConditionConfigBlockSize) {
			.block_size = 2,
		},
		// (HypergraphConfigFromFile) {
		// 	.file_format = MATRIX_MARKET,
		// 	.file_path = "1138_bus.mtx",
		// }
		(HypergraphConfigManual) {
			.vertex_count = 8,
			.hyperedge_indices = example_hyperedge_indices,
			.hyperedges = example_hyperedges,
		}
	);

	auto result = ctx.run_multi_node();

	printf("Runhund Final Result (2 lines):\n");
	print_vector(result.row_permutation);
	print_vector(result.column_permutation);

	MPI_Finalize();
}
