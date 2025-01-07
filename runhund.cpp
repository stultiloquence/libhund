#include <cstdio>
#include <memory>
#include <iostream>
#include <vector>
#include <thread>

#include <mtkahypar.h>
#include "context.h"



void try_mt_kahypar() {

	// Initialize thread pool
	mt_kahypar_initialize(
		std::thread::hardware_concurrency() /* use all available cores */,
		true /* activate interleaved NUMA allocation policy */
	);

	// Setup partitioning context
	mt_kahypar_context_t* context = mt_kahypar_context_new();
	mt_kahypar_load_preset(context, DEFAULT /* corresponds to MT-KaHyPar-D */);
	// In the following, we partition a hypergraph into two blocks
	// with an allowed imbalance of 3% and optimize the connective metric (KM1)
	mt_kahypar_set_partitioning_parameters(
		context,
		2 /* number of blocks */,
		0.03 /* imbalance parameter */,
		KM1 /* objective function */
	);
	mt_kahypar_set_seed(42 /* seed */);
	// Enable logging
	mt_kahypar_set_context_parameter(context, VERBOSE, "1");

	// Load Hypergraph for DEFAULT preset


	size_t example_hyperedge_indices[] = { 0, 3, 6, 8 };
	// long int example_hyperedges[] = { 0, 4, 7, 2, 4, 7, 0, 7 };
	long int example_hyperedges[] = { 0, 2, 3, 1, 2, 3, 0, 3 };
	// size_t example_hyperedge_indices[] = { 0, 3, 5, 8, 11, 15, 17, 19, 22, 24 };
	// long int example_hyperedges[] = { 0, 4, 7, 1, 3, 2, 4, 7, 3, 5, 6, 1, 2, 3, 4, 1, 5, 5, 6, 2, 3, 7, 0, 7 };
	mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(DEFAULT, 4, 3, example_hyperedge_indices, (mt_kahypar_hyperedge_id_t*) example_hyperedges, nullptr, nullptr);

	// Partition Hypergraph
	mt_kahypar_partitioned_hypergraph_t partitioned_hg = mt_kahypar_partition(hypergraph, context);

	// Extract Partition
	std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
	mt_kahypar_get_partition(partitioned_hg, partition.get());

	// Extract Block Weights
	std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
	mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

	// Compute Metrics
	const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
	const double km1 = mt_kahypar_km1(partitioned_hg);

	// Output Results
	std::cout << "Partition:" << std::endl;
	for (int i = 0; i < 4; i++) {
		std::cout << partition.get()[i] << ", ";
	}
	std::cout << std::endl;
	std::cout << "Partitioning Results:" << std::endl;
	std::cout << "Imbalance         = " << imbalance << std::endl;
	std::cout << "Km1               = " << km1 << std::endl;
	std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
	std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;

	mt_kahypar_free_context(context);
	mt_kahypar_free_hypergraph(hypergraph);
	mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}

int main(int argc, char const *argv[]) {
	//try_mt_kahypar();
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
		(HypergraphConfigManual) {
			.vertex_count = 8,
			.hyperedge_indices = example_hyperedge_indices,
			.hyperedges = example_hyperedges,
		}
	);

	auto result = ctx.run();

	print_vector(result.row_permutation);
	print_vector(result.column_permutation);
}
