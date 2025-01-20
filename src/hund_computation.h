#pragma once

#include <cfloat>
#include <climits>
#include <variant>
#include <vector>
#include <stdio.h>

#include <mtkahypar.h>
#include <mtkahypartypes.h>
#include <mpi.h>

#include <hypergraph.h>
#include <types.h>
#include <hypergraph_bisection.h>
#include <util.h>
#include <logger.h>
#include <initialize_mt_kahypar.h>

struct Partition {
	std::vector<int> partition;
	double quality; // Despite the name, smaller qualities indicate a better bisection
	double true_imbalance;
};

class HUNDComputation {

private:
	// Configuration Objects
	BisectionConfig bisection_config;
	BreakConditionConfig break_condition_config;
	Logger &logger;

	mt_kahypar_context_t* mt_kahypar_context;

	MPI_Comm MPI_COMM_ROOT;
	int mpi_node_id;
	int mpi_nr_of_nodes;

	// Hypergraph Data
	Hypergraph &hypergraph;

	Partition bisect(
		Hypergraph hypergraph
	) {
		BisectionConfigMtKahypar *cmk = std::get_if<BisectionConfigMtKahypar>(&bisection_config);
		if (cmk) {
			auto mt_hypergraph = mt_kahypar_create_hypergraph(
				DETERMINISTIC,
				hypergraph.get_vertex_count(),
				hypergraph.get_edge_count(),
				hypergraph.hyperedge_indices.data(),
				hypergraph.hyperedges.data(),
				nullptr,
				nullptr
			);
			auto partitioned_hg = mt_kahypar_partition(mt_hypergraph, mt_kahypar_context);

			auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(mt_hypergraph));
			mt_kahypar_get_partition(partitioned_hg, partition.get());

			auto result = std::vector(partition.get(), partition.get() + hypergraph.get_vertex_count());
			double quality;
			switch (cmk->objective_function) {
			case KahyparObjectiveFunction::KM1:
				quality = mt_kahypar_km1(partitioned_hg);
				break;
			case KahyparObjectiveFunction::CUT:
				quality = mt_kahypar_cut(partitioned_hg);
				break;
			case KahyparObjectiveFunction::SOED:
				quality = mt_kahypar_soed(partitioned_hg);
				break;
			default:
				assert(false);
			}

			auto true_imbalance = mt_kahypar_imbalance(partitioned_hg, mt_kahypar_context);
			mt_kahypar_free_hypergraph(mt_hypergraph);
			mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
			return {
				.partition = result,
				.quality = quality,
				.true_imbalance = true_imbalance,
			};
		}

		assert(false); // Should not happen;
		return {};
	}

	bool break_condition_reached(
		Hypergraph hypergraph,
		int recursion_depth
	) {
		if (std::holds_alternative<BreakConditionConfigRecursionDepth>(break_condition_config)) {
			auto bcc = std::get<BreakConditionConfigRecursionDepth>(break_condition_config);
			return recursion_depth >= bcc.depth;
		} else if (std::holds_alternative<BreakConditionConfigBlockSize>(break_condition_config)) {
			auto bcc = std::get<BreakConditionConfigBlockSize>(break_condition_config);
			return hypergraph.get_vertex_count() <= bcc.max_block_size_inclusive || hypergraph.get_edge_count() <= bcc.max_block_size_inclusive;
		} else {
			// should not happen
			assert(false);
			return true;
		}
	}

	RowColPermutation run_single_node(
		Hypergraph hypergraph,
		int recursion_depth
	) {
		if (break_condition_reached(hypergraph, recursion_depth)) {
			return {
				.row_permutation = identity_permutation(hypergraph.get_vertex_count()),
				.column_permutation = identity_permutation(hypergraph.get_edge_count()),
			};
		}

		auto partition = HUNDComputation::bisect(hypergraph);

		logger.log_potential_partition(mpi_node_id, hypergraph, recursion_depth, mpi_node_id, mpi_node_id, bisection_config, partition.partition, partition.quality, partition.true_imbalance);

		auto hb = HypergraphBisection(partition.partition, hypergraph);

		logger.log_best_bisection(mpi_node_id, hypergraph, recursion_depth, mpi_node_id, mpi_node_id, hb);

		auto result_0 = run_single_node(hb.get_hypergraph_0(), recursion_depth + 1);
		auto result_1 = run_single_node(hb.get_hypergraph_1(), recursion_depth + 1);

		return {
			.row_permutation = compose_permutations(
				hb.get_vertex_permutation(), 
				combine_permutations(
					result_0.row_permutation,
					result_1.row_permutation
				)
			),
			.column_permutation = compose_permutations(
				hb.get_edge_permutation(),
				combine_permutations(
					result_0.column_permutation,
					result_1.column_permutation,
					identity_permutation(hb.get_separator_size())
				)
			)
		};
	}

	RowColPermutation run_multi_node(
		Hypergraph hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		MPI_Comm RANGE_COMM
	) {
		if (break_condition_reached(hypergraph, recursion_depth)) {
			// todo: compute local reordering using some library
			return {
				.row_permutation = identity_permutation(hypergraph.get_vertex_count()),
				.column_permutation = identity_permutation(hypergraph.get_edge_count()),
			};
		}

		if (range_end - range_start == 1) {
			return run_single_node(hypergraph, recursion_depth);
		}

		// Compute different and find best partition on all nodes within RANGE_COMM. We use RANGE_COMM for communication of partitions using MPI_Allgather and MPI_Bcast. (To communicate the resulting permutations later, we use point-to-point communication in MPI_COMM_ROOT).

		// todo: 1. figure out your node's bisection parameter variation
		// 2. compute partition and quality
		auto partition = HUNDComputation::bisect(hypergraph);

		logger.log_potential_partition(mpi_node_id, hypergraph, recursion_depth, range_start, range_end, bisection_config, partition.partition, partition.quality, partition.true_imbalance);

		// 3. communicate qualities to all other nodes working on this bisection.
		int range_size = range_end - range_start;
		std::vector<double> qualities(range_size);
		MPI_Allgather(&partition.quality, 1, MPI_DOUBLE,
			&qualities[0], 1, MPI_DOUBLE, RANGE_COMM);
		// 4. Figure out the node with the best partition (deterministically!). Smaller qualities indicate a better bisection
		double best_quality = DBL_MAX;
		size_t best_node_id_within_RANGE_COMM = 0;
		for (size_t i = 0; i < range_size; i++) {
			if (qualities[i] < best_quality) {
				best_quality = qualities[i];
				best_node_id_within_RANGE_COMM = i;
			}
		}
		// 5. The winning node sends its partition, the others receive. Overwrites partition.partition. Relies on the fact that all partitions are of the same size.
		MPI_Bcast(&partition.partition[0], partition.partition.size(), MPI_INT, best_node_id_within_RANGE_COMM, RANGE_COMM);

		auto hb = HypergraphBisection(partition.partition, hypergraph);

		logger.log_best_bisection(mpi_node_id, hypergraph, recursion_depth, range_start, range_end, hb);

		Hypergraph own_hypergraph, other_hypergraph;
		int sub_range_start, sub_range_end;

		int range_midpoint = range_end - (range_end - range_start) / 2;
		int lower_half_size = range_midpoint - range_start;
		bool is_lower_half = mpi_node_id < range_midpoint;
		if (is_lower_half) {
			sub_range_start = range_start;
			sub_range_end = range_midpoint;
			own_hypergraph = hb.get_hypergraph_0();
			other_hypergraph = hb.get_hypergraph_1();
		} else {
			sub_range_start = range_midpoint;
			sub_range_end = range_end;
			own_hypergraph = hb.get_hypergraph_1();
			other_hypergraph = hb.get_hypergraph_0();
		}

		MPI_Comm SUB_COMM;
		MPI_Comm_split(RANGE_COMM, is_lower_half, mpi_node_id, &SUB_COMM);
		auto sub_result = run_multi_node(
			own_hypergraph,
			recursion_depth + 1,
			sub_range_start,
			sub_range_end,
			SUB_COMM
		);

		// Concatenate permutations to create a single payload for MPI
		auto own_permutations = sub_result.row_permutation;
		own_permutations.insert(
			own_permutations.end(),
			sub_result.column_permutation.begin(),
			sub_result.column_permutation.end()
		);

		// Reserve space for other permutations
		int other_permutations_size = other_hypergraph.get_edge_count() + other_hypergraph.get_vertex_count();
		std::vector<unsigned long> other_permutations(other_permutations_size);
		
		// printf("mpi_node_id = %i,\nrecursion_depth = %i,\nrange_start = %i,\nrange_end = %i,\nsub_range_start = %i,\nsub_range_end = %i,\nown_hypergraph.get_vertex_count() = %li,\nown_hypergraph.get_edge_count() = %li,\nown_permutations.size() = %li,\nother_hypergraph.get_vertex_count() = %li,\nother_hypergraph.get_edge_count() = %li,\nother_permutations.size() = %li\n", mpi_node_id, recursion_depth, range_start, range_end, sub_range_start, sub_range_end, own_hypergraph.get_vertex_count(), own_hypergraph.get_edge_count(), own_permutations.size(), other_hypergraph.get_vertex_count(), other_hypergraph.get_edge_count(), other_permutations.size());

		bool is_uneven = (range_end - range_start) % 2 == 1;
		bool is_last_in_lower_half = mpi_node_id + 1 == range_midpoint;
		bool is_last_in_upper_half = mpi_node_id + 1 == range_end;

		if (!(is_uneven && is_last_in_lower_half)) {
			int partner = is_lower_half ? mpi_node_id + lower_half_size : mpi_node_id - lower_half_size;

			MPI_Request send_request;

			// printf("[Node %d] talking to %d, sending %zu entries, ready to receive %zu entries.\n", mpi_node_id, partner, own_permutations.size(), other_permutations.size());
			// printf("[Node %d] Some addresses: %p %p %p", mpi_node_id, (void*) &own_permutations[0], (void*) &own_permutations, (void*) &send_request);

			MPI_Isend(&own_permutations[0], own_permutations.size(), MPI_UNSIGNED_LONG,
				partner, 1, MPI_COMM_ROOT, &send_request);
			MPI_Recv(&other_permutations[0], other_permutations_size, MPI_UNSIGNED_LONG,
				partner, 1, MPI_COMM_ROOT, MPI_STATUS_IGNORE);
			MPI_Wait(&send_request, MPI_STATUS_IGNORE);
		}

		if (is_uneven && is_last_in_upper_half) {
			MPI_Send(&own_permutations[0], own_permutations.size(), MPI_UNSIGNED_LONG,
				range_midpoint - 1, 1, MPI_COMM_ROOT);
		}

		if (is_uneven && is_last_in_lower_half) {
			MPI_Recv(&other_permutations[0], other_permutations_size, MPI_UNSIGNED_LONG,
				range_end - 1, 1, MPI_COMM_ROOT, MPI_STATUS_IGNORE);
		}

		// printf("mpi_node_id = %i has own permutation and other permutation:\n", mpi_node_id);
		// print_vector(own_permutations);
		// print_vector(other_permutations);

		std::vector<unsigned long> row_permutation, column_permutation;
		if (is_lower_half) {
			row_permutation = combine_permutations(
				&own_permutations[0],
				own_hypergraph.get_vertex_count(),
				&other_permutations[0],
				other_hypergraph.get_vertex_count()
			);
			column_permutation = combine_permutations(
				&own_permutations[own_hypergraph.get_vertex_count()],
				own_hypergraph.get_edge_count(),
				&other_permutations[other_hypergraph.get_vertex_count()],
				other_hypergraph.get_edge_count(),
				&identity_permutation(hb.get_separator_size())[0],
				hb.get_separator_size()
			);
		} else {
			row_permutation = combine_permutations(
				&other_permutations[0],
				other_hypergraph.get_vertex_count(),
				&own_permutations[0],
				own_hypergraph.get_vertex_count()
			);
			column_permutation = combine_permutations(
				&other_permutations[other_hypergraph.get_vertex_count()],
				other_hypergraph.get_edge_count(),
				&own_permutations[own_hypergraph.get_vertex_count()],
				own_hypergraph.get_edge_count(),
				&identity_permutation(hb.get_separator_size())[0],
				hb.get_separator_size()
			);
		}

		return {
			.row_permutation = compose_permutations(
				hb.get_vertex_permutation(),
				row_permutation
			),
			.column_permutation = compose_permutations(
				hb.get_edge_permutation(),
				column_permutation
			)
		};
	}


public:
	HUNDComputation(
		BisectionConfig bc,
		BreakConditionConfig bcc,
		MultithreadingConfig mc,
		Hypergraph &hypergraph,
		Logger &logger
	) : logger(logger), hypergraph(hypergraph) {
		this->bisection_config = bc;
		this->break_condition_config = bcc;

		BisectionConfigMtKahypar *cmk = std::get_if<BisectionConfigMtKahypar>(&bc);
		if (cmk) {
			initialize_mt_kahypar(mc.number_of_threads_per_rank);
			mt_kahypar_context = mt_kahypar_context_new();
			mt_kahypar_load_preset(mt_kahypar_context, DEFAULT);
			mt_kahypar_set_partitioning_parameters(
				mt_kahypar_context,
				2 /* number of blocks */,
				cmk->max_imbalance,
				to_mt_kahypar_objective_function(cmk->objective_function)
			);
			mt_kahypar_set_context_parameter(mt_kahypar_context, VERBOSE, "0");
		}

		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_ROOT);
	    MPI_Comm_size(MPI_COMM_ROOT, &mpi_nr_of_nodes);
	    MPI_Comm_rank(MPI_COMM_ROOT, &mpi_node_id);

	    assert(1 <= mpi_nr_of_nodes);
	    assert(0 <= mpi_node_id && mpi_node_id < mpi_nr_of_nodes);
	}

	~HUNDComputation() {
		if (std::holds_alternative<BisectionConfigMtKahypar>(bisection_config)) {
			mt_kahypar_free_context(mt_kahypar_context);
		}
	};

	RowColPermutation run_single_node() {
		return run_single_node(hypergraph, 0);
	}

	RowColPermutation run_multi_node() {
		return run_multi_node(hypergraph, 0, 0, mpi_nr_of_nodes, MPI_COMM_ROOT);
	}

	Hypergraph getHypergraph() {
		return hypergraph;
	}
};