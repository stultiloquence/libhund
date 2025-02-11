#pragma once
 
#include <mtkahypar.h>
#include <mtkahypartypes.h>
#include <hypergraph.h>
#include <initialize_mt_kahypar.h>

/**
 * This class represents a bisection of a hypergraph by Mt-KaHyPar into
 * possibly more than 2 blocks. The size of the separator can be computed to be
 * compared to other dissections. For more details on Mt-KaHyPar see [their
 * documentation](https://github.com/kahypar/mt-kahypar).
 */
class KahyparComputation {
private:
	mt_kahypar_context_t* mt_kahypar_context;
	Hypergraph &hypergraph;
	double cut;
	double km1;
	double soed;
public:
	/**
	 * Constructs an instance of an Mt-KaHyPar bisection. This initializes 
	 * everything required for the actual bisection. For more details on the 
	 * setup and parameters see [their documentation](https://github.com/kahypar/mt-kahypar).

	 * @param max_imbalance Maximum imbalance of the two partitions of the
	 *     bisection.
	 * @param number_of_threads Number of threads used for multi-threading. A
	 *     value of 0 uses the maximum number available.
	 * @param nr_of_blocks Number of blocks to dissect the hypergraph into.
	 * @param hypergraph The hypergraph to be dissected. 
	 */
	KahyparComputation(
		double max_imbalance,
		int number_of_threads,
		int nr_of_blocks,
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

	/**
	 * Destructor that makes sure that the Mt-KaHyPar library exits 
	 * cleanly.
	 */
	~KahyparComputation() {
		mt_kahypar_free_context(mt_kahypar_context);
	};

	/**
	 * Computes a dissection of the Hypergraph with Mt-KaHyPar with the 
	 * parameters set in the constructor. Also computes the three different 
	 * metrics that are offered by the library.
	 *
	 * @return std::vector<int> Partition of the vertices. This is a list of
	 * the same length as the hypergraph has vertices. The value for each index
	 * indicates which partition that vertex has been dissected into.
	 */
	std::vector<int> dissect() {
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
		cut = mt_kahypar_cut(partitioned_hg);
		km1 = mt_kahypar_km1(partitioned_hg);
		soed = mt_kahypar_soed(partitioned_hg);

		auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(mt_hypergraph));
		mt_kahypar_get_partition(partitioned_hg, partition.get());

		auto result = std::vector(partition.get(), partition.get() + hypergraph.get_vertex_count());

		mt_kahypar_free_hypergraph(mt_hypergraph);
		mt_kahypar_free_partitioned_hypergraph(partitioned_hg);

		return result;
	}

	/**
	 * @return double Get the value of the Connectivity Metric (KM1).
	 */
	double get_km1() {
		return km1;
	}
	
	/**
	 * @return double Get the value of the Cut-Net Metric.
	 */
	double get_cut() {
		return cut;
	}

	/**
	 * @return double Get the value of the Sum-of-external-degrees Metric.
	 */
	double get_soed() {
		return soed;
	}
};