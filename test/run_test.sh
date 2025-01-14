#!/bin/bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../external/mt-kahypar/build/lib
export LD_LIBRARY_PATH
mpiexec -np 2 ./test