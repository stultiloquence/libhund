#pragma once

#include <mtkahypar.h>
#include <thread>

bool __mt_kahypar_initialized = false;

void initialize_mt_kahypar() {
	if (!__mt_kahypar_initialized) {
		mt_kahypar_initialize(
			std::thread::hardware_concurrency() /* use all available cores */,
			true /* activate interleaved NUMA allocation policy */
		);
		__mt_kahypar_initialized = true;
		
		mt_kahypar_set_seed(42);
	}
}