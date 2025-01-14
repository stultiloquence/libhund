#pragma once

#include "types.h"
#include <cassert>
#include <iostream>
#include <vector>

template <typename T>
void print_vector(std::vector<T> v) {
	for (size_t i = 0; i < v.size(); i++) {
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

std::vector<unsigned long> identity_permutation(
	unsigned long size
) {
	std::vector<unsigned long> result(size);
	for (unsigned long i = 0; i < size; i++) {
		result[i] = i;
	}
	return result;
}

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