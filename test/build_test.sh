#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
mpicxx -std=c++17 -I../external -I../external/mt-kahypar/include -I../src/ -pthread test.cpp ../external/mt-kahypar/build/lib/libmtkahypar.so ../external/catch2/catch_amalgamated.a -o test