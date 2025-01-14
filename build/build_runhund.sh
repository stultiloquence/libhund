#!/bin/bash
mpicxx -std=c++17 ../src/runhund.cpp ../external/mt-kahypar/build/lib/libmtkahypar.so -I../external/mt-kahypar/include -I../src/ -I../external -pthread -o runhund