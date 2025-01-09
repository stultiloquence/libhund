#pragma once

#include <cfloat>
#include <climits>
#include <filesystem>
#include <variant>
#include <vector>
#include <thread>
#include <fstream>
#include <stdio.h>

#include <mtkahypar.h>
#include <mtkahypartypes.h>
#include "fast_matrix_market/fast_matrix_market.hpp"
#include <mpi.h>

#include "types.h"
#include "hypergraph_bisection.h"
#include "util.h"

struct BisectionConfigMtKahypar {
	double max_imbalance;
	kahypar_objective_function_t objective_function;
};

typedef std::variant<BisectionConfigMtKahypar> BisectionConfig;

struct BreakConditionConfigRecursionDepth {
	int depth;
};

struct BreakConditionConfigBlockSize {
	int block_size;
};

typedef std::variant<BreakConditionConfigRecursionDepth, BreakConditionConfigBlockSize> BreakConditionConfig;

struct HypergraphConfigManual {
	size_t vertex_count;
	std::vector<size_t> hyperedge_indices;
	std::vector<unsigned long> hyperedges;
};

struct HypergraphConfigFromFile {
	matrix_file_format_t file_format;
	std::filesystem::path file_path;
};

typedef std::variant<HypergraphConfigManual, HypergraphConfigFromFile> HypergraphConfig;

struct Partition {
	std::vector<int> partition;
	double quality;
};

class Context {

private:
	// Configuration Objects
	BisectionConfig bisection_config;
	BreakConditionConfig break_condition_config;
	HypergraphConfig hypergraph_config;

	inline static bool mt_kahypar_initialized = false;
	mt_kahypar_context_t* mt_kahypar_context;

	inline static bool mpi_initialized = false;
	MPI_Comm MPI_COMM_ROOT;
	int mpi_node_id;
	int mpi_nr_of_nodes;

	// Hypergraph Data
	Hypergraph hypergraph;

	void initialize_mt_kahypar() {
		mt_kahypar_initialize(
			std::thread::hardware_concurrency() /* use all available cores */,
			true /* activate interleaved NUMA allocation policy */
		);
		mt_kahypar_initialized = true;
	}

	static mt_kahypar_objective_t to_mt_kahypar_objective_function(
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

	void sort_by_col_first(
		std::vector<unsigned long> rows,
		std::vector<unsigned long> cols,
		std::vector<unsigned long> &sorted_rows,
		std::vector<unsigned long> &sorted_cols
	) {
		// Find sort permutation
		auto perm = identity_permutation(rows.size());
		std::sort(perm.begin(), perm.end(),
			[&](std::size_t i, std::size_t j) {
				if (cols[i] != cols[j]) {
					return cols[i] < cols[j];
				}
				if (rows[i] != rows[j]) {
					return rows[i] < rows[j];
				}
				return false;
		 });

		// Apply permutation
		sorted_rows.reserve(rows.size());
		std::transform(perm.begin(), perm.end(), std::back_inserter(sorted_rows), [&](auto i) { return rows[i]; });
		sorted_cols.reserve(cols.size());
		std::transform(perm.begin(), perm.end(), std::back_inserter(sorted_cols), [&](auto i) { return cols[i]; });
	}

	int load_hypergraph_from_matrix_file(
		HypergraphConfigFromFile hcff
	) {
		if (hcff.file_format == MATRIX_MARKET) {			
			std::vector<unsigned long> rows, cols;
			std::vector<double> vals;
			unsigned long row_count = 0, col_count = 0;

			std::ifstream matrix_file(hcff.file_path);
			fast_matrix_market::read_matrix_market_triplet(
				matrix_file,
				row_count, col_count,
				rows, cols, vals);

			std::vector<unsigned long> sorted_rows, sorted_cols;
			sort_by_col_first(rows, cols, sorted_rows, sorted_cols);
	
			hypergraph.vertex_count = row_count;
			hypergraph.hyperedges = sorted_rows;

			hypergraph.hyperedge_indices.reserve(col_count + 1);
			hypergraph.hyperedge_indices.push_back(0);
			unsigned int i = 0;
			for (unsigned long col = 0; col < col_count; col++) {
				while (i < sorted_cols.size() && sorted_cols[i] <= col) {
					i++;
				}
				hypergraph.hyperedge_indices.push_back(i);
			}
			return 0;
		} else {
			throw std::runtime_error("Only MATRIX_MARKET format is supported right now.");
		}
	}

	int load_hypergraph_directly(
		HypergraphConfigManual hcm
	) {
		hypergraph = {
			.vertex_count = hcm.vertex_count,
			.hyperedges = hcm.hyperedges,
			.hyperedge_indices = hcm.hyperedge_indices,
		};
		return 0;
	}

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
			case HUND_KM1:
				quality = mt_kahypar_km1(partitioned_hg);
				break;
			case HUND_CUT:
				quality = mt_kahypar_cut(partitioned_hg);
				break;
			case HUND_SOED:
				quality = mt_kahypar_soed(partitioned_hg);
				break;
			default:
				assert(false);
			}

			mt_kahypar_free_hypergraph(mt_hypergraph);
			mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
			return {
				.partition = result,
				.quality = quality
			};
		}

		assert(false); // Should not happen;
		return {
			.partition = std::vector<int>(),
			.quality = 0.0,
		};
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
			return hypergraph.get_vertex_count() <= bcc.block_size || hypergraph.get_edge_count() <= bcc.block_size;
		} else {
			// should not happen
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

		auto partition = Context::bisect(hypergraph);
		auto hb = HypergraphBisection(partition.partition, hypergraph);

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
		auto partition = Context::bisect(hypergraph);
		// 3. communicate qualities to all other nodes working on this bisection.
		int range_size = range_end - range_start;
		std::vector<double> qualities(range_size);
		MPI_Allgather(&partition.quality, 1, MPI_DOUBLE,
			&qualities[0], 1, MPI_DOUBLE, RANGE_COMM);
		// 4. Figure out the node with the best partition (deterministically!)
		double best_quality = -DBL_MAX;
		size_t best_node_id_within_RANGE_COMM = 0;
		for (size_t i = 0; i < range_size; i++) {
			if (qualities[i] > best_quality) {
				best_quality = qualities[i];
				best_node_id_within_RANGE_COMM = i;
			}
		}
		// 5. The winning node sends its partition, the others receive. Overwrites partition.partition. Relies on the fact that all partitions are of the same size.
		MPI_Bcast(&partition.partition[0], partition.partition.size(), MPI_UNSIGNED_LONG, best_node_id_within_RANGE_COMM, RANGE_COMM);

		auto hb = HypergraphBisection(partition.partition, hypergraph);

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
		
		printf("mpi_node_id = %i,\nrecursion_depth = %i,\nrange_start = %i,\nrange_end = %i,\nsub_range_start = %i,\nsub_range_end = %i,\nown_hypergraph.get_vertex_count() = %li,\nown_hypergraph.get_edge_count() = %li,\nown_permutations.size() = %li,\nother_hypergraph.get_vertex_count() = %li,\nother_hypergraph.get_edge_count() = %li,\nother_permutations.size() = %li\n", mpi_node_id, recursion_depth, range_start, range_end, sub_range_start, sub_range_end, own_hypergraph.get_vertex_count(), own_hypergraph.get_edge_count(), own_permutations.size(), other_hypergraph.get_vertex_count(), other_hypergraph.get_edge_count(), other_permutations.size());

		bool is_uneven = (range_end - range_start) % 2 == 1;
		bool is_last_in_lower_half = mpi_node_id + 1 == range_midpoint;
		bool is_last_in_upper_half = mpi_node_id + 1 == range_end;

		if (!(is_uneven && is_last_in_lower_half)) {
			int partner = is_lower_half ? mpi_node_id + lower_half_size : mpi_node_id - lower_half_size;

			MPI_Request *send_request;
			MPI_Isend(&own_permutations[0], own_permutations.size(), MPI_UNSIGNED_LONG,
				partner, 1, MPI_COMM_ROOT, send_request);
			MPI_Recv(&other_permutations[0], other_permutations_size, MPI_UNSIGNED_LONG,
				partner, 1, MPI_COMM_ROOT, MPI_STATUS_IGNORE);
			MPI_Wait(send_request, MPI_STATUS_IGNORE);	
		}

		if (is_uneven && is_last_in_upper_half) {
			MPI_Send(&own_permutations[0], own_permutations.size(), MPI_UNSIGNED_LONG,
				range_midpoint - 1, 1, MPI_COMM_ROOT);
		}

		if (is_uneven && is_last_in_lower_half) {
			MPI_Recv(&other_permutations[0], other_permutations_size, MPI_UNSIGNED_LONG,
				range_end - 1, 1, MPI_COMM_ROOT, MPI_STATUS_IGNORE);
		}

		printf("mpi_node_id = %i has own permutation and other permutation:\n", mpi_node_id);
		print_vector(own_permutations);
		print_vector(other_permutations);

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
	Context(
		BisectionConfig bc,
		BreakConditionConfig bcc,
		HypergraphConfig hc
	) {
		this->bisection_config = bc;
		this->break_condition_config = bcc;
		this->hypergraph_config = hc;

		BisectionConfigMtKahypar *cmk = std::get_if<BisectionConfigMtKahypar>(&bc);
		if (cmk) {
			if (!mt_kahypar_initialized) {
				initialize_mt_kahypar();
			}
			mt_kahypar_set_seed(42);
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

		if (std::holds_alternative<HypergraphConfigManual>(hc)) {
			load_hypergraph_directly(std::get<HypergraphConfigManual>(hc));
		} else if (std::holds_alternative<HypergraphConfigFromFile>(hc)) {
			load_hypergraph_from_matrix_file(std::get<HypergraphConfigFromFile>(hc));
		}

		if (!mpi_initialized) {
			MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_ROOT);
		    MPI_Comm_size(MPI_COMM_ROOT, &mpi_nr_of_nodes);
		    MPI_Comm_rank(MPI_COMM_ROOT, &mpi_node_id);

		    assert(1 <= mpi_nr_of_nodes);
		    assert(0 <= mpi_node_id && mpi_node_id < mpi_nr_of_nodes);

			mpi_initialized = true;
		}

	}
	~Context() {
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