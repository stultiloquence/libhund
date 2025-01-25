#pragma once

#include <cassert>
#include <iostream>
#include <vector>

/**
 * Helper function to print a std::vector<T> to a std::ostream.
 */
template <typename T>
void print_vector_to_stream(std::vector<T> v, std::ostream &os) {
	for (size_t i = 0; i < v.size() - 1; i++) {
		os << v[i] << " ";
	}
	os << v[v.size() - 1] << std::endl;
}

/**
 * Helper function to print a std::vector<T> to std::cout.
 */
template <typename T>
void print_vector(std::vector<T> v) {
	print_vector_to_stream(v, std::cout);
}

/**
 * Construct a new identity permutation on n elements.
 * @param n The size of the permutation, that is, the number of elements
 *     being permuted.
 * @return a std::vector of length n containing, in order, the unsigned
 *     longs 0 through to n-1. Represents the identity permutation on
 *     n elements.
 */
std::vector<unsigned long> identity_permutation(
	unsigned long n
) {
	std::vector<unsigned long> result(n);
	for (unsigned long i = 0; i < n; i++) {
		result[i] = i;
	}
	return result;
}

/**
 * Combine two or three permutations into one permutation that applies
 * the first permutation onto the first first_size many elements, the
 * second permutation onto the next second_size many elements, and then
 * potentially the third on the next third_size many elements.
 * @param first, first_size Pointer to the first element and length of a
 *     C-style array containing the first permutation.
 * @param second, second_size Pointer to the second element and length of a
 *     C-style array containing the second permutation.
 * @param third, third_size Pointer to the third element and length of a
 *     C-style array containing the third permutation. Optional, defaults to
 *     an empty array on 0 elements.
 * @return a std::vector of length first_size + second_size + third_size
 *     representing a permutation that applies the first permutation onto
 *     the first first_size many elements, the second permutation onto the
 *     next second_size many elements, and then potentially the third on the
 *     next third_size many elements.
 */
std::vector<unsigned long> combine_permutations(
	unsigned long *first,
	size_t first_size,
	unsigned long *second,
	size_t second_size,
	unsigned long *third = nullptr,
	size_t third_size = 0
) {
	std::vector<unsigned long> result;
	result.reserve(first_size + second_size + third_size);
	for (size_t i = 0; i < first_size; i++) {
		result.push_back(first[i]);
	}
	for (size_t i = 0; i < second_size; i++) {
		result.push_back(second[i] + first_size);
	}
	for (size_t i = 0; i < third_size; i++) {
		result.push_back(third[i] + first_size + second_size);
	}
	return result;
}

/**
 * Combine two or three permutations into one permutation that applies
 * the first permutation onto the first first_size many elements, the
 * second permutation onto the next second_size many elements, and then
 * potentially the third on the next third_size many elements.
 * @param first vector containing the first permutation.
 * @param second vector containing the second permutation.
 * @param third vector containing the third permutation. Optional, defaults to
 *     an empty vector.
 * @return a std::vector of length first.size() + second.size() + third.size()
 *     representing a permutation that applies the first permutation onto
 *     the first first.size() many elements, the second permutation onto the
 *     next second.size() many elements, and then potentially the third on the
 *     next third.size() many elements.
 */
std::vector<unsigned long> combine_permutations(
	std::vector<unsigned long> first,
	std::vector<unsigned long> second,
	std::vector<unsigned long> third = std::vector<unsigned long>()
) {
	return combine_permutations(
		&first[0],
		first.size(),
		&second[0],
		second.size(),
		third.empty() ? nullptr : &third[0],
		third.size()
	);
}

/**
 * Computes the mathematical composition of two permutations (interpreted as
 * bijective functions on the set { 0, ..., n-1 }).
 * @param first The first permutation to be applied.
 * @param second The second permutation to be applied.
 * @return the permutation resulting from first applying the first and then
 *     the second permutation to the identity permutation. Equivalently, the
 *     composition second o first.
 */
std::vector<unsigned long> compose_permutations(
	std::vector<unsigned long> first,
	std::vector<unsigned long> second
) {
	assert(first.size() == second.size());
	std::vector<unsigned long> result;
	unsigned long length = first.size();
	result.reserve(length);
	for (unsigned long v = 0; v < length; v++) {
		result.push_back(second[first[v]]);
	}
	return result;
}