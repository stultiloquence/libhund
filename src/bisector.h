#pragma once

#include <vector>

#include <mtkahypar.h>
#include <mtkahypartypes.h>

#include <initialize_mt_kahypar.h>
#include <hypergraph.h>

/**
 * Represents a partition of a set of partition.size() vertices into two sets.
 * A Partition just holds a vector containing a partition id for each index.
 * For instance, partition.partition[3] == 1 signifies that the vertex with
 * id 3 belongs to partition 1. The struct also holds information about the
 * quality and imbalance of the partition.
 */
struct Partition {
	std::vector<int> partition;
	double quality; // Despite the name, smaller qualities indicate a better bisection
	double true_imbalance;
};

/**
 * Abstract base class for bisectors. Represents a single, fully configured
 * hypergraph bisection algorithm, that can be run repeatedly on different
 * hypergraphs.
 */
class Bisector {
public:
	/**
	 * Bisects a hypergraph, returning a vertex partition.
	 * @param hypergraph the hypergraph to be partitioned.
	 * @return A Partition instance indicating which partition (0 or 1) each
	 *     hypergraph vertex belongs to, along with quality and imbalance of
	 *     the partition.
	 */
	virtual Partition bisect(Hypergraph &hypergraph) = 0;
};

/**
 * Bisects hypergraphs using the MtKaHyPar library.
 */
class MtKahyparBisector : public Bisector {
private:
	mt_kahypar_context_t* mt_kahypar_context;
public:
	/**
	 * Fully configures the MtKaHyPar library. Initializes the library if that
	 * has not happened yet (using initialize_mt_kahypar).
	 * @param max_imbalance The maximum imbalance parameter passed along to
	 *     MtKaHyPar.
	 * @param number_of_threads Number of threads to be used. Defaults to 0,
	 *     which means using all available threads.
	 */
	MtKahyparBisector(
  		double max_imbalance,
		int number_of_threads_per_rank = 0
	) {
		initialize_mt_kahypar(number_of_threads_per_rank);
		mt_kahypar_context = mt_kahypar_context_new();
		mt_kahypar_load_preset(mt_kahypar_context, DEFAULT);
		mt_kahypar_set_partitioning_parameters(
			mt_kahypar_context,
			2 /* number of blocks */,
			max_imbalance,
			KM1
		);
		mt_kahypar_set_context_parameter(mt_kahypar_context, VERBOSE, "0");
	}

	/**
	 * Destructor, freeing the associated MtKaHyPar context.
	 */
	~MtKahyparBisector() {
		mt_kahypar_free_context(mt_kahypar_context);
	}

	/**
	 * Bisects the given hypergraph using the MtKaHyPar library.
	 * @param hypergraph the hypergraph to be partitioned.
	 * @return A Partition instance indicating which partition (0 or 1) each
	 *     hypergraph vertex belongs to, along with quality and imbalance of
	 *     the partition.
	 */
	Partition bisect(
		Hypergraph &hypergraph
	) override {
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
		double quality = mt_kahypar_km1(partitioned_hg);

		auto true_imbalance = mt_kahypar_imbalance(partitioned_hg, mt_kahypar_context);
		mt_kahypar_free_hypergraph(mt_hypergraph);
		mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
		return {
			.partition = result,
			.quality = quality,
			.true_imbalance = true_imbalance,
		};
	}
};