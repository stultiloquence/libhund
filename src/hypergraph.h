#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>
#include <fstream>

#include <fast_matrix_market/fast_matrix_market.hpp>

#include <util.h>

/**
 * Types of matrix files from which a hypergraph can be constructed using the 
 * Hypergraph constrcutor.
 */
typedef enum MatrixFileFormat {
	MATRIX_MARKET
	// MATLAB,
	// RUTHERFORD_BOEING,
} matrix_file_format_t;

/**
 * This class represents an instance of a hypergraph. It is represents by the
 * Compressed Column Storage of the underlying matrix without any values.
 */
class Hypergraph {
private:

	/**
	 * Sorting the entries from the Matrix Market format increasingly by columns 
	 *  and then rows.
	 *
	 * @param rows Vector of the row entries.
	 * @param cols Vector of the column entries.
	 * @param sorted_rows Reference to the vector in which to store the sorted
	 *     rows.
	 * @param sorted_cols Reference to the vector in which to store the sorted
	 *     columns. 
	 */
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
public:
	size_t vertex_count;
	std::vector<unsigned long> hyperedges;
	std::vector<size_t> hyperedge_indices;

	/**
	 * Construct a default hyergraph with no vertices or hypereges.
	 */
	Hypergraph() {
		Hypergraph(0, std::vector<size_t>(), std::vector<unsigned long>());
	}

	
	/**
	 * Construct a hyergraph manually by passing all the information directly.
	 * This can be seen as passing a matrix in Compressed Column Storage 
	 * format.
	 *
	 * @param vertex_count Number of vertices of the hypergraph / Number of 
	 *     rows of the matrix.
	 * @param hyperedge_indices Vector that indicates at which index a new  
	 *     hyeredge begins in the hyperedges paramater / Vector that indicates 
	 *     at which index a new column begins in the following parameter.
	 * @param hyperedges Vector that contains the vertices each hyperedges 
	 *     spans, concatinated / Vector that contains the row index of each 
	 *     entry of the matrix, when iterating columnwise.
	 */
	Hypergraph(
		size_t vertex_count,
		std::vector<size_t> hyperedge_indices,
		std::vector<unsigned long> hyperedges
	) : vertex_count(vertex_count),
		hyperedges(hyperedges),
		hyperedge_indices(hyperedge_indices) {};

	/**
	 * Constructs a hypergraph from a (sparse) matrix in a given file format.
	 *
	 * @param file_format Format of the matrix file, e.g. Matrix Market.
	 * @param file_path Path to the file containg the matrix.
	 */
	Hypergraph(
		matrix_file_format_t file_format,
		std::filesystem::path file_path
	) {
		if (file_format == MATRIX_MARKET) {      
			std::vector<unsigned long> rows, cols;
			std::vector<double> vals;
			unsigned long row_count = 0, col_count = 0;

			std::ifstream matrix_file(file_path);
			fast_matrix_market::read_matrix_market_triplet(
				matrix_file,
				row_count, col_count,
				rows, cols, vals);

			std::vector<unsigned long> sorted_rows, sorted_cols;
			sort_by_col_first(rows, cols, sorted_rows, sorted_cols);
	
			this->vertex_count = row_count;
			this->hyperedges = sorted_rows;

			this->hyperedge_indices.reserve(col_count + 1);
			this->hyperedge_indices.push_back(0);
			unsigned int i = 0;
			for (unsigned long col = 0; col < col_count; col++) {
				while (i < sorted_cols.size() && sorted_cols[i] <= col) {
					i++;
				}
				this->hyperedge_indices.push_back(i);
			}
		} else {
			// No other format is supported as of now.
    		assert(false);
		}
	}

	/**
	 * @return size_t Get the number of vertices in the hyerpgraph. 
	 */
	size_t get_vertex_count() {
		return vertex_count;
	}

	/**
	 * @return size_t Get the number of hyperedges in the hyerpgraph. 
	 */
	size_t get_edge_count() {
		return hyperedge_indices.size() - 1;
	}

};