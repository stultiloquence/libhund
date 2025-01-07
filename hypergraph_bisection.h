#pragma once

#include "types.h"
#include <iostream>

class HypergraphBisection {

private:
	Hypergraph original_hypergraph;

	size_t vertex_count;
	size_t edge_count;

	std::vector<int> vertex_partition;
	unsigned long vertex_partition_sizes[2] = { 0, 0 };
	std::vector<unsigned long> vertex_position_within_partition;
	
	// Partition id that each edge belongs to. A value of 2 indicates that the edge belongs to neither partition 0 or 1 completely.
	std::vector<int> edge_partition;
	unsigned long edge_partition_sizes[3] = { 0, 0, 0 };
	std::vector<unsigned long> edge_position_within_partition;

	Hypergraph hypergraph_0;
	Hypergraph hypergraph_1;
	unsigned long separator_size;

	std::vector<unsigned long> vertex_permutation;
	std::vector<unsigned long> edge_permutation;

	void compute_vertex_vectors() {
		vertex_position_within_partition.reserve(vertex_count);

		for (size_t i = 0; i < vertex_count; i++) {
			unsigned long next_position = vertex_partition_sizes[vertex_partition[i]];
			vertex_position_within_partition.push_back(next_position);
			vertex_partition_sizes[vertex_partition[i]] += 1;
		}
	}

	void compute_edge_vectors() {
		edge_partition.reserve(edge_count);
		edge_position_within_partition.reserve(edge_count);

		for (size_t i = 0; i < edge_count; i++) {
			size_t edge_start = original_hypergraph.hyperedge_indices[i];
			size_t edge_end = original_hypergraph.hyperedge_indices[i + 1];

			unsigned long first_vertex = original_hypergraph.hyperedges[edge_start];
			
			int hyperedge_partition = vertex_partition[first_vertex];
			
			for (size_t j = edge_start + 1; j < edge_end; j++) {
				unsigned long vertex = original_hypergraph.hyperedges[j];
				if (hyperedge_partition != vertex_partition[vertex]) {
					hyperedge_partition = 2;
					break;
				}
			}

			edge_partition.push_back(hyperedge_partition);
			edge_position_within_partition.push_back(edge_partition_sizes[hyperedge_partition]);
			edge_partition_sizes[hyperedge_partition] += 1;
		}
	}

	void compute_sub_hypergraphs() {
		for (size_t edge = 0; edge < edge_count; edge++) {
			size_t edge_start = original_hypergraph.hyperedge_indices[edge];
			size_t edge_end = original_hypergraph.hyperedge_indices[edge + 1];

			if (edge_partition[edge] == 2) {
				// Edges with edge_partition[edge] == 2 go into neither subgraph.
				continue;
			}

			Hypergraph *hg = (edge_partition[edge] == 0) ? &hypergraph_0 : &hypergraph_1;
			hg->hyperedge_indices.push_back(hg->hyperedges.size());

			for (size_t j = edge_start; j < edge_end; j++) {
				unsigned long vertex = original_hypergraph.hyperedges[j];
				hg->hyperedges.push_back(vertex_position_within_partition[vertex]);
			}
		}

		hypergraph_0.hyperedge_indices.push_back(hypergraph_0.hyperedges.size());
		hypergraph_1.hyperedge_indices.push_back(hypergraph_1.hyperedges.size());

		hypergraph_0.vertex_count = vertex_partition_sizes[0];
		hypergraph_1.vertex_count = vertex_partition_sizes[1];

		separator_size = edge_count
			- (hypergraph_0.get_edge_count())
			- (hypergraph_1.get_edge_count());
	}

	unsigned long single_vertex_absolute_position(
		unsigned long vertex
	) {
		int partition = vertex_partition[vertex];
		if (partition == 0) {
			return vertex_position_within_partition[vertex];
		} else {
			return vertex_partition_sizes[0] + vertex_position_within_partition[vertex];
		}
	}

	unsigned long single_edge_absolute_position(
		unsigned long edge
	) {
		int partition = edge_partition[edge];
		if (partition == 0) {
			return edge_position_within_partition[edge];
		} else if (partition == 1) {
			return edge_partition_sizes[0] + edge_position_within_partition[edge];
		} else {
			return edge_partition_sizes[0] + edge_partition_sizes[1] + edge_position_within_partition[edge];
		}
	}

	void compute_permutations() {
		vertex_permutation.reserve(vertex_count);
		for (unsigned long vertex = 0; vertex < vertex_count; vertex++) {
			vertex_permutation.push_back(single_vertex_absolute_position(vertex));
		}
		edge_permutation.reserve(edge_count);
		for (unsigned long edge = 0; edge < edge_count; edge++) {
			edge_permutation.push_back(single_edge_absolute_position(edge));
		}
	}

public:
	HypergraphBisection(
		std::vector<int> vert_partition,
		Hypergraph hypergraph
	) {
		original_hypergraph = hypergraph;

		vertex_count = original_hypergraph.get_vertex_count();
		edge_count = original_hypergraph.get_edge_count();
		vertex_partition = vert_partition;

		compute_vertex_vectors();
		compute_edge_vectors();
		compute_sub_hypergraphs();
		compute_permutations();
	};

	Hypergraph get_hypergraph_0() {
		return hypergraph_0;
	};

	Hypergraph get_hypergraph_1() {
		return hypergraph_1;
	};

	std::vector<unsigned long> get_vertex_permutation() {
		return vertex_permutation;
	}

	std::vector<unsigned long> get_edge_permutation() {
		return edge_permutation;
	}

	unsigned long get_separator_size() {
		return separator_size;
	}

};