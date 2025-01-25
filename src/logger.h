#pragma once

#include <algorithm>
#include <mpi.h>
#include <vector>

#include <mpi.h>

#include <types.h>
#include <bisector.h>
#include <hypergraph_bisection.h>
#include <hypergraph.h>

/**
 * Abstract class for Loggers that can be passed to the HUNDComputation class.
 */
class Logger {
public:
	/**
	 * This logging function is called in both HUNDComputation.run_single_node
	 * and .run_multi_node, right after a bisection is performed. It informs
	 * the logger over the computed bisection. It is a "potential" partition
	 * because in the multi node case, it might not be the best partition out
	 * of all the parallel nodes attempting the same bisection. So this logging
	 * function can in particular be used to compare different node's partition
	 * attempts.
	 * @param node_id The logging node's MPI rank.
	 * @param hypergraph The (sub)hypergraph that is currently being partitioned.
	 * @param recursion_depth The current recursion depth.
	 * @param range_start The first (inclusive) MPI rank of all the ranks
	 *     that are attempting to bisect the same hypergraph.
	 * @param range_end The first (exclusive) MPI rank of all the ranks
	 *     that are attempting to bisect the same hypergraph.
	 * @param partition The computed partition.
	 */
	virtual void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		Partition &partition
	) {}

	/**
	 * This logging function is called in both HUNDComputation.run_single_node
	 * and .run_multi_node, after all partition attempts have been compared,
	 * the best one has been chosen, and applied to the hypergraph, resulting in
	 * the HypergraphBisection hb. It only differes from log_potential_partition
	 * in the multi node case, and only when there is more than one node attempting
	 * the partition.
	 * @param node_id The logging node's MPI rank.
	 * @param hypergraph The (sub)hypergraph that is currently being partitioned.
	 * @param recursion_depth The current recursion depth.
	 * @param range_start The first (inclusive) MPI rank of all the ranks
	 *     that are attempting to bisect the same hypergraph.
	 * @param range_end The first (exclusive) MPI rank of all the ranks
	 *     that are attempting to bisect the same hypergraph.
	 * @param hb The computed hypergraph bisection.
	 */
	virtual void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) {}
};

/**
 * This logger does nothing. To be passed to HUNDComputation when no logging
 * is desired. Incurs no performance penalty.
 */
class NoOpLogger : public Logger {
public:
	/**
	 * Does nothing.
	 */
	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		Partition &partition
	) override {}
	/**
	 * Does nothing.
	 */
	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override {}
};

/**
 * SeparatorSizeLogger. Records the total separator size that the HUND
 * algorithm yields. Also computes a total separator size weighted by
 * the number of blocks each separator spans. Can only be used when the
 * HUNDComputation is configured with a BreakConditionConfigRecursionDepth.
 */
class SeparatorSizeLogger : public Logger {
private:
	int total_recursion_depth;
	unsigned long total_separator_size_weighted = 0;
	unsigned long total_separator_size = 0;
public:
	/**
	 * To compute the weighted partition, the Logger needs to know
	 * the maximum recursion depth.
	 */
	SeparatorSizeLogger(
		int total_recursion_depth
	) : total_recursion_depth(total_recursion_depth) {};

	/**
	 * Does nothing.
	 */
	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		Partition &partition
	) override {}

	/**
	 * Adds the next bisection's separator size to the ongoing count.
	 */
	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override {
		if (node_id == range_start) {
			total_separator_size_weighted += hb.get_separator_size() * std::pow(2, total_recursion_depth - recursion_depth);
			total_separator_size += hb.get_separator_size();
		}
	}

	/**
	 * Makes all MPI nodes collect each other's data, so that they all contain
	 * the correct final results. To be called after the HUNDComputation has
	 * finished.
	 */
	void gather() {
		MPI_Comm MPI_COMM_SEPARATOR_SIZE;
		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_SEPARATOR_SIZE);


		int world_size, node_id;
		MPI_Comm_size(MPI_COMM_SEPARATOR_SIZE, &world_size);
		MPI_Comm_rank(MPI_COMM_SEPARATOR_SIZE, &node_id);

		std::vector<unsigned long> total_sizes_weighted(world_size);

		MPI_Allgather(&total_separator_size_weighted, 1, MPI_UNSIGNED_LONG,
			&total_sizes_weighted[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_SEPARATOR_SIZE);
		total_separator_size_weighted = 0;
		for (auto size : total_sizes_weighted) {
			total_separator_size_weighted += size;
		}

		std::vector<unsigned long> total_sizes(world_size);
		MPI_Allgather(&total_separator_size, 1, MPI_UNSIGNED_LONG,
			&total_sizes[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_SEPARATOR_SIZE);
		total_separator_size = 0;
		for (auto size : total_sizes) {
			total_separator_size += size;
		}
	}

	/**
	 * Return the total weighted separator size. This is the sum of the number of
	 * columns spanned by each separator, multiplied by the number of blocks the
	 * separator spans. Currently uses the maximum number of blocks potentially
	 * spanned by each separator, and is therefore rather broken.
	 * @return The total weighted separator size.
	 */
	unsigned long get_total_separator_size_weighted() {
		return total_separator_size_weighted;
	}

	/**
	 * Return the total separator size. This is the sum of the number of
	 * columns spanned by each separator.
	 * @return The total separtor size.
	 */
	unsigned long get_total_separator_size() {
		return total_separator_size;
	}
};

/**
 * Records information about a single BisectionAttempt that occured during
 * the HUNDComputation algorithm.
 */
struct BisectionAttempt {
	int node_id;
	int range_start;
	int range_end;
	size_t vertex_count;
	size_t edge_count;
	double quality;
	double true_imbalance;
	size_t separator_size;
	bool operator==(const BisectionAttempt &rhs) {
		return this->node_id == rhs.node_id
			&& this->range_start == rhs.range_start
			&& this->range_end == rhs.range_end
			&& this->vertex_count == rhs.vertex_count
			&& this->edge_count == rhs.edge_count
			&& this->quality == rhs.quality
			&& this->true_imbalance == rhs.true_imbalance
			&& this->separator_size == rhs.separator_size;
	}
};

/**
 * Records aggregated information about all bisection attempts that were
 * done during the HUNDComputation for one given subhypergraph, that is, at
 * a given recursion depth by the nodes with node id [range_start, range_end[,
 * all attempting the same bisection. This can be used to compare the different
 * bisection attempts.
 */
struct BisectionAttempts {
	int recursion_depth;
	int range_start;
	int range_end;
	size_t vertex_count;
	size_t edge_count;
	std::vector<double> qualities;
	std::vector<double> true_imbalances;
	std::vector<size_t> separator_sizes;
};

/**
 * Writes a std::vector<double> to an std::ostream in JSON format.
 * @param v the vector to be written to the stream.
 * @param os the stream to be written into.
 */
void print_as_json(std::vector<double> v, std::ostream &os) {
	if (v.size() == 0) {
		os << "[]";
		return;
	}
	os << "[";
	for (size_t i = 0; i < v.size() - 1; i++) {
		os << v[i] << ", ";
	}
	os << v[v.size() - 1] << "]";
}

/**
 * Write a std::vector<size_t> to an std::ostream in JSON format.
 * @param v the vector to be written to the stream.
 * @param os the stream to be written into.
 */
void print_as_json(std::vector<size_t> v, std::ostream &os) {
	if (v.size() == 0) {
		os << "[]";
		return;
	}
	os << "[";
	for (size_t i = 0; i < v.size() - 1; i++) {
		os << v[i] << ", ";
	}
	os << v[v.size() - 1] << "]";
}

/**
 * Write a BisectionAttempts to an std::ostream in JSON format.
 * @param bisection_attempts The BisectionAttempts object to be serialized.
 * @param os the stream to be written into.
 */
void print_as_json(BisectionAttempts bisection_attempts, std::ostream &os) {
	os << "{\n";
	os << "    \"recursion_depth\": " << bisection_attempts.recursion_depth << ",\n";
	os << "    \"range_start\": " << bisection_attempts.range_start << ",\n";
	os << "    \"range_end\": " << bisection_attempts.range_end << ",\n";
	os << "    \"vertex_count\": " << bisection_attempts.vertex_count << ",\n";
	os << "    \"edge_count\": " << bisection_attempts.edge_count << ",\n";
	os << "    \"qualities\": ";
	print_as_json(bisection_attempts.qualities, os);
	os << ",\n";
	os << "    \"true_imbalances\": ";
	print_as_json(bisection_attempts.true_imbalances, os);
	os << ",\n";
	os << "    \"separator_sizes\": ";
	print_as_json(bisection_attempts.separator_sizes, os);
	os << "\n}";
}

/**
 * Write a std::vector<BisectionAttempts> to an std::ostream in JSON format.
 * @param bas The vector of BisectionAttempts to be serialized.
 * @param os the stream to be written into.
 */
void print_as_json(std::vector<BisectionAttempts> &bas, std::ostream &os) {
	if (bas.size() == 0) {
		os << "[]";
		return;
	}
	os << "[\n";
	for (size_t i = 0; i < bas.size() - 1; i++) {
		print_as_json(bas[i], os);
		os << ", ";
	}
	print_as_json(bas[bas.size() - 1], os);
	os << "\n]";
}

/**
 * Logger that collects all the attempts at parallel bisection that occured
 * during the HUNDComputation as a list. Can be used to see the effect that
 * the parallel bisection attempts have.
 */
class BisectionQualityRangeLogger : public Logger {
private:
	constexpr static const BisectionAttempt null_value = {
		-1, -1, -1, 0, 0, 0.0, 0.0, 0
	};
	int node_id = -1;
	std::vector<BisectionAttempt> bisection_attempts_by_recursion_depth;
	MPI_Datatype MPI_BISECTION_ATTEMPT_TYPE;
public:
	/**
	 * Constructor. Prepares MPI stuff for logging.
	 */
	BisectionQualityRangeLogger() {
		// Create a type for struct BisectionAttempt.
		const int nitems = 8;
		int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1, 1};
		MPI_Datatype types[nitems] = {
			MPI_INT,
			MPI_INT,
			MPI_INT,
			MPI_UNSIGNED_LONG,
			MPI_UNSIGNED_LONG,
			MPI_DOUBLE,
			MPI_DOUBLE,
			MPI_UNSIGNED_LONG
		};
		MPI_Aint offsets[nitems];

		offsets[0] = offsetof(BisectionAttempt, node_id);
		offsets[1] = offsetof(BisectionAttempt, range_start);
		offsets[2] = offsetof(BisectionAttempt, range_end);
		offsets[3] = offsetof(BisectionAttempt, vertex_count);
		offsets[4] = offsetof(BisectionAttempt, edge_count);
		offsets[5] = offsetof(BisectionAttempt, quality);
		offsets[6] = offsetof(BisectionAttempt, true_imbalance);
		offsets[7] = offsetof(BisectionAttempt, separator_size);

		MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_BISECTION_ATTEMPT_TYPE);
		MPI_Type_commit(&MPI_BISECTION_ATTEMPT_TYPE);
	}

	/**
	 * Stores a potential partition as a BisectionAttempt if there are at least
	 * two nodes attempting this bisection in parallel.
	 */
	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		Partition &partition
	) override {
		if (range_start == range_end) {
			// There is only one bisection attempt.
			return;
		}

		if (this->node_id == -1) {
			this->node_id = node_id;
		} else if (this->node_id != node_id) {
			// There is one logger instance per node, so each instance should be called by one node_id only.
			assert(false);
		}
		auto hb = HypergraphBisection(partition.partition, hypergraph);
		bisection_attempts_by_recursion_depth.resize(std::max(
			(unsigned long) recursion_depth + 1,
			bisection_attempts_by_recursion_depth.size()
		));
		bisection_attempts_by_recursion_depth[recursion_depth] = {
			.node_id = node_id,
			.range_start = range_start,
			.range_end = range_end,
			.vertex_count = hypergraph.get_vertex_count(),
			.edge_count = hypergraph.get_edge_count(),
			.quality = partition.quality,
			.true_imbalance = partition.true_imbalance,
			.separator_size = hb.get_separator_size(),
		};
	}

	/**
	 * Does nothing.
	 */
	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override { }

	/**
	 * Gathers the information of all BisectionQualityRangeLoggers working on
	 * the different MPI nodes into one vector of BisectionAttempts, one for
	 * each bisection attempt that occured with at least two nodes attempting
	 * a bisection.
	 * @return A std::vector of BisectionAttempts.
	 */
	std::vector<BisectionAttempts> gather() {
		MPI_Comm MPI_COMM;
		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM);
		int world_size, rank;
		MPI_Comm_size(MPI_COMM, &world_size);
		MPI_Comm_rank(MPI_COMM, &rank);

		// 1. Find maximum recursion depth across all nodes
		size_t own_recursion_depth = bisection_attempts_by_recursion_depth.size();
		std::vector<size_t> recursion_depths = std::vector<size_t>(world_size);
		MPI_Allgather(&own_recursion_depth, 1, MPI_UNSIGNED_LONG, &recursion_depths[0], 1, MPI_UNSIGNED_LONG, MPI_COMM);
		size_t max_recursion_depth = *std::max_element(recursion_depths.begin(), recursion_depths.end());

		// 2. Pad own array with dummy values.
		for (auto i = bisection_attempts_by_recursion_depth.size(); i < max_recursion_depth; i++) {
			bisection_attempts_by_recursion_depth.push_back(null_value);
		}

		// 3. Gather each node's bisection_attempts_by_recursion_depths array.
		auto all_attempts = std::vector<BisectionAttempt>(world_size * max_recursion_depth);
		
		MPI_Allgather(&bisection_attempts_by_recursion_depth[0], max_recursion_depth, MPI_BISECTION_ATTEMPT_TYPE, &all_attempts[0], max_recursion_depth, MPI_BISECTION_ATTEMPT_TYPE, MPI_COMM);

		// 4. Transform the data
		auto attempts2d = std::vector<std::vector<BisectionAttempt>>(world_size);
		for (int i = 0; i < world_size; i++) {
			attempts2d[i] = std::vector<BisectionAttempt>(
				all_attempts.begin() + i * max_recursion_depth,
				all_attempts.begin() + (i + 1) * max_recursion_depth
			);
		}
		std::sort(attempts2d.begin(), attempts2d.end(), [](std::vector<BisectionAttempt> &a, std::vector<BisectionAttempt> &b) {
			return a[0].node_id < b[0].node_id;
		});

		auto result = std::vector<BisectionAttempts>();
		auto next = BisectionAttempts{};
		for (auto depth = 0; depth < max_recursion_depth; depth++) {
			for (auto id = 0; id < world_size; id++) {
				auto entry = attempts2d[id][depth];
				if (entry == null_value) {
					continue;
				}
				next.qualities.push_back(entry.quality);
				next.true_imbalances.push_back(entry.true_imbalance);
				next.separator_sizes.push_back(entry.separator_size);
				if (entry.range_end == entry.node_id + 1) {
					next.recursion_depth = depth;
					next.range_start = entry.range_start;
					next.range_end = entry.range_end;
					next.vertex_count = entry.vertex_count;
					next.edge_count = entry.edge_count;
					result.push_back(next);
					next = BisectionAttempts{};
				}
			}
		}
		return result;
	}
};