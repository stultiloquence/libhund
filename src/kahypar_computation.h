
#include <mtkahypar.h>
#include <mtkahypartypes.h>
#include <types.h>
#include <hypergraph.h>
#include <initialize_mt_kahypar.h>

class KahyparComputation {
private:
	mt_kahypar_context_t* mt_kahypar_context;
	Hypergraph &hypergraph;
	double km1;
public:
	KahyparComputation(
		double max_imbalance,
		int nr_of_blocks,
		int number_of_threads,
		Hypergraph &hypergraph
	) : hypergraph(hypergraph) {
		initialize_mt_kahypar(number_of_threads);

		mt_kahypar_context = mt_kahypar_context_new();
		mt_kahypar_load_preset(mt_kahypar_context, DEFAULT);
		mt_kahypar_set_partitioning_parameters(
			mt_kahypar_context,
			nr_of_blocks,
			max_imbalance,
			KM1
		);
		mt_kahypar_set_context_parameter(mt_kahypar_context, VERBOSE, "0");
	}

	std::vector<int> bisect() {
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
		km1 = mt_kahypar_km1(partitioned_hg);

		auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(mt_hypergraph));
		mt_kahypar_get_partition(partitioned_hg, partition.get());

		auto result = std::vector(partition.get(), partition.get() + hypergraph.get_vertex_count());

		mt_kahypar_free_hypergraph(mt_hypergraph);
		mt_kahypar_free_partitioned_hypergraph(partitioned_hg);

		return result;
	}

	RowColPermutation run(std::vector<int>) {
		// todo
		return {};
	}

	RowColPermutation run() {
		return run(bisect());
	}

	unsigned long size_of_separator(
		std::vector<int> partition
	) {
		unsigned long result = 0;
		auto edge_count = hypergraph.get_edge_count();
		for (size_t i = 0; i < edge_count; i++) {
			size_t edge_start = hypergraph.hyperedge_indices[i];
			size_t edge_end = hypergraph.hyperedge_indices[i + 1];

			if (edge_start == edge_end) {
				continue;
			}

			unsigned long first_vertex = hypergraph.hyperedges[edge_start];
			int hyperedge_partition = partition[first_vertex];
			
			for (size_t j = edge_start + 1; j < edge_end; j++) {
				unsigned long vertex = hypergraph.hyperedges[j];
				if (hyperedge_partition != partition[vertex]) {
					result += 1;
					break;
				}
			}
		}
		return result;
	}

	unsigned long size_of_separator() {
		return size_of_separator(bisect());
	}

	double get_km1() {
		return km1;
	}

	~KahyparComputation() {
		mt_kahypar_free_context(mt_kahypar_context);
	};
};