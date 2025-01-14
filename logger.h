#pragma once

#include <mpi.h>
#include <vector>

#include <mpi.h>

#include <types.h>
#include <hypergraph_bisection.h>
#include <hypergraph.h>

class Logger {
public:
	virtual void log_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		std::vector<int> &partition,
		double partition_quality
	) {}
	virtual void log_bisection(
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
	void log_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		std::vector<int> &partition,
		double partition_quality
	) override {}
	void log_bisection(
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
	void log_partition(
		int node_id,
		Hypergraph &hypergraph,
		int recursion_depth,
		int range_start,
		int range_end,
		std::vector<int> &partition,
		double partition_quality
	) override {}

	void log_bisection(
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