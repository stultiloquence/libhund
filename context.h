#pragma once

#include <climits>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <variant>
#include <vector>
#include <thread>
#include <fstream>

#include <mtkahypar.h>
#include <mtkahypartypes.h>
#include "fast_matrix_market/fast_matrix_market.hpp"

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

class Context {

private:
	// Configuration Objects
	BisectionConfig bisection_config;
	BreakConditionConfig break_condition_config;
	HypergraphConfig hypergraph_config;

	inline static bool mt_kahypar_initialized = false;
	mt_kahypar_context_t* mt_kahypar_context;

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

	std::vector<int> bisect(
		Hypergraph hypergraph
	) {
		if (std::holds_alternative<BisectionConfigMtKahypar>(bisection_config)) {

			auto mt_hypergraph = mt_kahypar_create_hypergraph(
				DEFAULT,
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

			mt_kahypar_free_hypergraph(mt_hypergraph);
			mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
			return result;
		}

		assert(false); // Should not happen;
		return std::vector<int>();
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

	RowColPermutation run(
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
		auto hb = HypergraphBisection(partition, hypergraph);

		auto result_0 = run(hb.get_hypergraph_0(), recursion_depth + 1);
		auto result_1 = run(hb.get_hypergraph_1(), recursion_depth + 1);

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

	}
	~Context() {
		if (std::holds_alternative<BisectionConfigMtKahypar>(bisection_config)) {
			mt_kahypar_free_context(mt_kahypar_context);
		}
	};

	RowColPermutation run() {
		return run(hypergraph, 0);
	}

	Hypergraph getHypergraph() {
		return hypergraph;
	}
};