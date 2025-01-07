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
	std::vector<unsigned long> first,
	std::vector<unsigned long> second,
	std::vector<unsigned long> third = std::vector<unsigned long>()
) {
	std::vector<unsigned long> result;
	result.reserve(first.size() + second.size() + third.size());
	for (auto v : first) {
		result.push_back(v);
	}
	for (auto v : second) {
		result.push_back(v + first.size());
	}
	for (auto v : third) {
		result.push_back(v + first.size() + second.size());
	}
	return result;
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