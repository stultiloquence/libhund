#pragma once

#include <mtkahypar.h>
#include <thread>


bool __mt_kahypar_initialized = false;

/**
 * Helper function to initilaize Mt-KaHyPar to make sure the library doesn't 
 * get instatiated multiple times, e.g. during the separator size test.
 */
void initialize_mt_kahypar(
	int nr_of_threads
) {
	if (!__mt_kahypar_initialized) {
		mt_kahypar_initialize(
			(nr_of_threads == 0) ? std::thread::hardware_concurrency() : nr_of_threads,
			true /* activate interleaved NUMA allocation policy */
		);
		__mt_kahypar_initialized = true;
		
		mt_kahypar_set_seed(42);
	}
}