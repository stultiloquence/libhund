#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>
#include <fstream>

#include <fast_matrix_market/fast_matrix_market.hpp>

#include <util.h>

typedef enum MatrixFileFormat {
	// MATLAB,
	// RUTHERFORD_BOEING,
	MATRIX_MARKET
} matrix_file_format_t;

class Hypergraph {
private:
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

	size_t get_vertex_count() {
		return vertex_count;
	}

	size_t get_edge_count() {
		return hyperedge_indices.size() - 1;
	}

	Hypergraph() {
		Hypergraph(0, std::vector<size_t>(), std::vector<unsigned long>());
	}

	Hypergraph(
		size_t vertex_count,
		std::vector<size_t> hyperedge_indices,
		std::vector<unsigned long> hyperedges
	) : vertex_count(vertex_count),
		hyperedges(hyperedges),
		hyperedge_indices(hyperedge_indices) {};

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
			throw std::runtime_error("Only MATRIX_MARKET format is supported right now.");
		}
	}

};