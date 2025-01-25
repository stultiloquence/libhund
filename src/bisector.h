#include <vector>

#include <mtkahypar.h>
#include <mtkahypartypes.h>

#include <initialize_mt_kahypar.h>
#include <hypergraph.h>


struct Partition {
	std::vector<int> partition;
	double quality; // Despite the name, smaller qualities indicate a better bisection
	double true_imbalance;
};

class Bisector {
public:
	virtual Partition bisect(Hypergraph hypergraph) = 0;
};

class MtKahyparBisector : public Bisector {
private:
	mt_kahypar_context_t* mt_kahypar_context;
public:
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

	~MtKahyparBisector() {
		mt_kahypar_free_context(mt_kahypar_context);
	}

	Partition bisect(
		Hypergraph hypergraph
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