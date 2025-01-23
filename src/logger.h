#pragma once

#include <algorithm>
#include <mpi.h>
#include <vector>

#include <mpi.h>

#include <types.h>
#include <hypergraph_bisection.h>
#include <hypergraph.h>

class Logger {
public:
	virtual void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		BisectionConfig bisection_config,
		std::vector<int> &partition,
		double partition_quality,
		double true_imbalance
	) {}
	virtual void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) {}
};

class NoOpLogger : public Logger {
public:
	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		BisectionConfig bisection_config,
		std::vector<int> &partition,
		double partition_quality,
		double true_imbalance
	) override {}
	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override {}
};

class SeparatorSizeLogger : public Logger {
private:
	unsigned long total_separator_size = 0;
public:
	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		BisectionConfig bisection_config,
		std::vector<int> &partition,
		double partition_quality,
		double true_imbalance
	) override {}

	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override {
		if (node_id == range_start) {
			total_separator_size += hb.get_separator_size();
		}
	}

	void gather() {
		MPI_Comm MPI_COMM_SEPARATOR_SIZE;
		MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_SEPARATOR_SIZE);


		int world_size, node_id;
		MPI_Comm_size(MPI_COMM_SEPARATOR_SIZE, &world_size);
		MPI_Comm_rank(MPI_COMM_SEPARATOR_SIZE, &node_id);

		std::vector<unsigned long> total_sizes(world_size);

		MPI_Allgather(&total_separator_size, 1, MPI_UNSIGNED_LONG,
			&total_sizes[0], 1, MPI_UNSIGNED_LONG, MPI_COMM_SEPARATOR_SIZE);

		total_separator_size = 0;
		for (auto size : total_sizes) {
			total_separator_size += size;
		}
	}

	unsigned long get_total_separator_size() {
		return total_separator_size;
	}
};

struct BisectionAttempt {
	int node_id;
	int range_start;
	int range_end;
	double max_imbalance;
	double quality;
	double true_imbalance;
	double relative_separator_size;
	bool operator==(const BisectionAttempt &rhs) {
		return this->node_id == rhs.node_id
			&& this->range_start == rhs.range_start
			&& this->range_end == rhs.range_end
			&& this->max_imbalance == rhs.max_imbalance
			&& this->quality == rhs.quality
			&& this->true_imbalance == rhs.true_imbalance
			&& this->relative_separator_size == rhs.relative_separator_size;
	}
};

struct BisectionAttempts {
	int recursion_depth;
	int range_start;
	int range_end;
	std::vector<double> max_imbalances;
	std::vector<double> qualities;
	std::vector<double> true_imbalances;
	std::vector<double> relative_separator_sizes;
};

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

void print_as_json(BisectionAttempts bisection_attempts, std::ostream &os) {
	os << "{\n";
	os << "    \"recursion_depth\": " << bisection_attempts.recursion_depth << ",\n";
	os << "    \"range_start\": " << bisection_attempts.range_start << ",\n";
	os << "    \"range_end\": " << bisection_attempts.range_end << ",\n";
	os << "    \"max_imbalances\": ";
	print_as_json(bisection_attempts.max_imbalances, os);
	os << ",\n";
	os << "    \"qualities\": ";
	print_as_json(bisection_attempts.qualities, os);
	os << ",\n";
	os << "    \"true_imbalances\": ";
	print_as_json(bisection_attempts.true_imbalances, os);
	os << ",\n";
	os << "    \"relative_separator_sizes\": ";
	print_as_json(bisection_attempts.relative_separator_sizes, os);
	os << "\n}";
}

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

class BisectionQualityRangeLogger : public Logger {
private:
	constexpr static const BisectionAttempt null_value = {
		-1, -1, -1, 0.0, 0.0, 0.0, 0.0
	};
	int node_id = -1;
	std::vector<BisectionAttempt> bisection_attempts_by_recursion_depth;
	MPI_Datatype MPI_BISECTION_ATTEMPT_TYPE;
public:
	BisectionQualityRangeLogger() {
		// Create a type for struct BisectionAttempt.
		const int nitems = 7;
		int blocklengths[nitems] = {1, 1, 1, 1, 1, 1, 1};
		MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
		MPI_Aint offsets[nitems];

		offsets[0] = offsetof(BisectionAttempt, node_id);
		offsets[1] = offsetof(BisectionAttempt, range_start);
		offsets[2] = offsetof(BisectionAttempt, range_end);
		offsets[3] = offsetof(BisectionAttempt, max_imbalance);
		offsets[4] = offsetof(BisectionAttempt, quality);
		offsets[5] = offsetof(BisectionAttempt, true_imbalance);
		offsets[6] = offsetof(BisectionAttempt, relative_separator_size);

		MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_BISECTION_ATTEMPT_TYPE);
		MPI_Type_commit(&MPI_BISECTION_ATTEMPT_TYPE);
	}

	void log_potential_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		BisectionConfig bisection_config,
		std::vector<int> &partition,
		double partition_quality,
		double true_imbalance
	) override {
		BisectionConfigMtKahypar *cmk = std::get_if<BisectionConfigMtKahypar>(&bisection_config);
		if (cmk) {
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
			auto hb = HypergraphBisection(partition, hypergraph);
			bisection_attempts_by_recursion_depth.resize(std::max(
				(unsigned long) recursion_depth + 1,
				bisection_attempts_by_recursion_depth.size()
			));
			bisection_attempts_by_recursion_depth[recursion_depth] = {
				.node_id = node_id,
				.range_start = range_start,
				.range_end = range_end,
				.max_imbalance = cmk->max_imbalance,
				.quality = partition_quality,
				.true_imbalance = true_imbalance,
				.relative_separator_size = ((double) hb.get_separator_size()) / hypergraph.get_edge_count(),
			};
		}
	}

	void log_best_bisection(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		HypergraphBisection &hb
	) override { }

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
				next.max_imbalances.push_back(entry.max_imbalance);
				next.qualities.push_back(entry.quality);
				next.true_imbalances.push_back(entry.true_imbalance);
				next.relative_separator_sizes.push_back(entry.relative_separator_size);
				if (entry.range_end == entry.node_id + 1) {
					next.recursion_depth = depth;
					next.range_start = entry.range_start;
					next.range_end = entry.range_end;
					result.push_back(next);
					next = BisectionAttempts{};
				}
			}
		}
		return result;
	}
};