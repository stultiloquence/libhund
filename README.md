# libhund

## Installation

libhund can be used either locally or on the HPC cluster.

### Locally

### HPC Cluster

```
module load GCC/11.3.0
module load Boost/1.79.0
module load tbb/2021.5.0
module load hwloc/2.7.1
module load OpenMPI/4.1.4

# install https://github.com/kahypar/mt-kahypar into mt-kahypar. build as a library. see their readme.md

git clone --depth=1 https://github.com/alugowski/fast_matrix_market

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:mt-kahypar/build/lib $MPICXX -std=c++17 -DNDEBUG runhund.cpp -o runhund mt-kahypar/build/lib/libmtkahypar.so -Imt-kahypar/include -I. -pthread

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:mt-kahypar/build/lib $MPIEXEC $FLAGS_MPI_BATCH ./runhund
```

### Testing

(what we did once:)
```
1. download catch_amalgamated.cpp and .hpp into catch2/
2. patch catch_amalgamated.cpp according to https://stackoverflow.com/a/58290117
3. mpicxx -std=c++17 -DCATCH_AMALGAMATED_CUSTOM_MAIN -c catch2/catch_amalgamated.cpp
4. ar rvs catch_amalgamated.a catch_amalgamated.o
5. 

```

every time:

```
mpicxx -std=c++17 test.cpp catch_amalgamated.a -o test -I.
mpiexec -np 4 ./test // ignore random output from multiple nodes

mpiexec -np 1 ./test [1rank]
mpiexec -np 2 ./test [2rank]
```