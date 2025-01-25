#!/usr/bin/zsh

#SBATCH --time=00:02:00 # max runtime

### Number of nodes.
#SBATCH --nodes=2
#SBATCH --partition=devel
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=3
#SBATCH --job-name=hundcli_fd18
#SBATCH --output=hundcli_fd18.txt

### Beginning of executable commands
module load GCC/11.3.0 Boost/1.79.0 tbb/2021.5.0 hwloc/2.7.1 OpenMPI/4.1.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../external/mt-kahypar/build/lib

$MPIEXEC -n 32 ../build/hundcli --threads-per-rank=3 \
    --matrix-file=../examples/fd18.mtx --matrix-file-format MATRIX_MARKET \
    --bisection-method=MT_KAHYPAR --kahypar-max-imbalance=0.05 \
    --break-condition=RECURSION_DEPTH --recursion-depth=3 \
    run --output-type=FILES --row-file=rows.txt --col-file=cols.txt