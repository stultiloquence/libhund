#pragma once

#include <hypergraph.h>


/** 
 * This class represents a single bisection of a hypergraph. It computes the 
 * new sub-hypergraphs and corresponding permutations based on a hypergraph
 * and the partition of its vertices.
 */
class HypergraphBisection {

private:
	static const int partition_for_empty_edge = 0;

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

	/**
	 * Fills the two vectors that state how big each vertex partition is and 
	 * what position each vertex has within its own partition.   
	 */
	void compute_vertex_vectors() {
		vertex_position_within_partition.reserve(vertex_count);

		for (size_t i = 0; i < vertex_count; i++) {
			unsigned long next_position = vertex_partition_sizes[vertex_partition[i]];
			vertex_position_within_partition.push_back(next_position);
			vertex_partition_sizes[vertex_partition[i]] += 1;
		}
	}

	/**
	 * Fills three vectors with information about the hyperedges. They state
	 * which partition each edge belongs to or if the edge is in the separator,
	 * how big each hyperedge partition is and what position each vertex has
	 * within its own partition. 
	 *
	 * For this it has to check whether each hyperedge is comepletly contained
	 * in one partition or belongs to the separator.
	 */
	void compute_edge_vectors() {
		edge_partition = std::vector<int>(edge_count);
		edge_position_within_partition =  std::vector<unsigned long>(edge_count);

		for (size_t i = 0; i < edge_count; i++) {
			size_t edge_start = original_hypergraph.hyperedge_indices[i];
			size_t edge_end = original_hypergraph.hyperedge_indices[i + 1];

			if (edge_start == edge_end) {
				// Empty hyperedges don't make much sense in the context of a hypergraph. However, our hypergraphs arise from matrices, where empty columns might exist. It is unclear where a empty hyperedge belongs, so we arbitrarily choose partition_for_empty_edge.
				edge_partition[i] = partition_for_empty_edge;
				edge_position_within_partition[i] = edge_partition_sizes[partition_for_empty_edge];
				edge_partition_sizes[partition_for_empty_edge] += 1;
				continue;
			}

			unsigned long first_vertex = original_hypergraph.hyperedges[edge_start];
			
			int hyperedge_partition = vertex_partition[first_vertex];
			
			for (size_t j = edge_start + 1; j < edge_end; j++) {
				unsigned long vertex = original_hypergraph.hyperedges[j];
				if (hyperedge_partition != vertex_partition[vertex]) {
					hyperedge_partition = 2;
					break;
				}
			}

			edge_partition[i] = hyperedge_partition;
			edge_position_within_partition[i] = edge_partition_sizes[hyperedge_partition];
			edge_partition_sizes[hyperedge_partition] += 1;
		}
	}

	/**
	 * Creates the the two sub-hypergraphs based on the vertex and hyperedge
	 * vectors and determines the size of the separator with this information.
	 *
	 * Checks the partition of each edge and creates the new hyeredges for the
	 * appropriate sub-hypergraph with the information where each vertex is
	 * relative within this sub-hypergraph.   
	 */
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

	/**
	 * Computes the absolute position in the original hypergraph for a vertex 
	 * based on the partition it is in and its size within this partition. 

	 * @param vertex Vertex whose position gets determined.  
	 * @return unsigned long Absolute position of the vertex.
	 */
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

	/**
	 * Computes the absolute position in the original hypergraph for a edge 
	 * based on the partition it is in and its size within this partition. 

	 * @param edge Hyperedge whose position gets determined.  
	 * @return unsigned long Absolute position of the hyperedge.
	 */
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

	/**
	 * Compute the permutation vectors that turn the matrix, which hyergraph is
	 * based upon, into the block structure computed from the bisection.
	 */
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
	/**
	 * Constructs an bisection of a hypergraph based on a permutation of the
	 * vertices. Computes all values that can possibly be queried in the 
	 * getters.  
	 *
	 * @param vert_partition Describes the partition of the vertices. The size
	 *     has to be the same as the amound of vertices of the hypergraph and
	 *     the values should only be 0 or 1. 
	 * @param hypergraph The hypergraph to be bisected.
	 */
	HypergraphBisection(
		std::vector<int> vert_partition,
		Hypergraph hypergraph
	) : original_hypergraph(hypergraph) {
		vertex_count = original_hypergraph.get_vertex_count();
		edge_count = original_hypergraph.get_edge_count();
		vertex_partition = vert_partition;

		compute_vertex_vectors();
		compute_edge_vectors();
		compute_sub_hypergraphs();
		compute_permutations();
	};

	/**
	 * @return Get the first sub-hypergraph created from the bisection. 
	 */
	Hypergraph get_hypergraph_0() {
		return hypergraph_0;
	};

	
	/**
	 * @return Get the second sub-hypergraph created from the bisection. 
	 */
	Hypergraph get_hypergraph_1() {
		return hypergraph_1;
	};

	
	/**
	 * @return Get the permuation vector for the rows of the matrix, which the
	 * vertices of the hypergraph are based upon, computed from the bisection.
	 */
	std::vector<unsigned long> get_vertex_permutation() {
		return vertex_permutation;
	}

	/**
	 * @return Get the permuation vector for the columns of the matrix, which 
	 * the hyperedges of the hypergraph are based upon, computed from the
	 * bisection.
	 */
	std::vector<unsigned long> get_edge_permutation() {
		return edge_permutation;
	}

	/**
	 * @return Get the number of edges that are part of the separator. 
	 */
	unsigned long get_separator_size() {
		return separator_size;
	}

};