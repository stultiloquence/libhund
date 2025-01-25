#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
mpicxx -std=c++17 ../src/hund_cli.cpp ../external/mt-kahypar/build/lib/libmtkahypar.so -I../external/mt-kahypar/include -I../src/ -I../external -pthread -o hundcli