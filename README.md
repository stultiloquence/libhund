# libhund

## Installation

libhund can be used either locally or on the HPC cluster.

### Locally

### HPC Cluster

```
module load GCC/11.3.0 Boost/1.79.0 tbb/2021.5.0 hwloc/2.7.1 OpenMPI/4.1.4

# install https://github.com/kahypar/mt-kahypar into mt-kahypar. build as a library. see their readme.md

git clone --depth=1 https://github.com/alugowski/fast_matrix_market

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:mt-kahypar/build/lib $MPICXX -std=c++17 -DNDEBUG runhund.cpp -o runhund mt-kahypar/build/lib/libmtkahypar.so -Imt-kahypar/include -I. -pthread

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:mt-kahypar/build/lib $MPIEXEC $FLAGS_MPI_BATCH ./runhund
```

### Testing

(what we did once:)
```
1. download catch_amalgamated.cpp and .hpp into catch2/
2. mpicxx -std=c++17 -DCATCH_AMALGAMATED_CUSTOM_MAIN -c catch2/catch_amalgamated.cpp
3. ar rvs catch_amalgamated.a catch_amalgamated.o
```

every time:

```
mpirun -n 2 ./run_test.sh [correctness]
```
