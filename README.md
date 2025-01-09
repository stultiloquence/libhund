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