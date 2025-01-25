# libhund

## Installation

libhund can be used either locally or on the HPC cluster.

First, you need to install all dependencies. libhund needs to be built from source as we do not distribute any binaries.

1. You need to install the following libraries with the given or a compatible versions.
```
GCC/11.3.0
Boost/1.79.0
tbb/2021.5.0
hwloc/2.7.1
OpenMPI/4.1.4
```
On the RWTH HPC compute cluster, this is accomplished by executing the command `module load GCC/11.3.0 Boost/1.79.0 tbb/2021.5.0 hwloc/2.7.1 OpenMPI/4.1.4`.

2. Clone the repository: `git clone https://stultiloquence/libhund && cd libhund`.

3. Install [mt-kahypar](https://github.com/kahypar/mt-kahypar) library into `external/mt-kahypar` following their instructions. Briefly:
	1. `cd external`
	2. `git clone https://github.com/kahypar/mt-kahypar`
	3. `cd mt-kahypar && mkdir build && cd build`
	4. `cmake .. --preset=default`
	5. `make install-mtkahypar`
	6. `cd ../../..` (return to the libhund root directory)

Now you can use the library through including the files in the `src/` directory.

## CLI Usage

### Building the CLI

After completing the steps under Installation, you can build the CLI through running the script under `build/build_hund_cli.sh`. This will place an executable `hundcli` into the `build/` folder.

### Using the CLI Locally

Since the `mt-kahypar` library is dynamically linked, you first need to tell the linker where the library can be found:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(realpath ./external/mt-kahypar/build/lib/)
```

Then, you can run the cli by executing `build/hundcli`. Here is an example command:

```
mkdir results
build/hundcli \
	--threads-per-rank=4 \
    --matrix-file=examples/paper_example.mtx --matrix-file-format=MATRIX_MARKET \
    --bisection-method=MT_KAHYPAR --kahypar-max-imbalance=0.05 \
    --break-condition=BLOCK_SIZE --block-size=2 \
    run --output-type=FILES --row-file=results/rows.txt --col-file=results/cols.txt
```

To run `hundcli` on multiple MPI ranks, use `mpiexec`:

```
mpiexec -n 4 build/hundcli \
	--threads-per-rank=4 \
    --matrix-file=examples/paper_example.mtx --matrix-file-format=MATRIX_MARKET \
    --bisection-method=MT_KAHYPAR --kahypar-max-imbalance=0.05 \
    --break-condition=BLOCK_SIZE --block-size=2 \
    run --output-type=FILES --row-file=results/rows.txt --col-file=results/cols.txt
```

In some situations (like when passing `--help`), every MPI rank outputs to stdout. To suppress this, use the `build/hund_cli_one_output.sh` helper script:

```
mpiexec -n 4 build/hund_cli_one_output.sh \
	--help
```

### Run on RWTH HPC Compute Cluster

To run `hundcli` on the RWTH HPC compute cluster, you can either run it like you would locally on the login node, in which case you first need to run the following two commands once per session:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(realpath ./external/mt-kahypar/build/lib/)
module load GCC/11.3.0 Boost/1.79.0 tbb/2021.5.0 hwloc/2.7.1 OpenMPI/4.1.4
```

and then invoke `hundcli` like you would locally.

To instead run the algorithm as a SLURM job, you can use and modify the example SLURM script `deploy/example_submit.sh`.

### CLI Options

todo: list of all options / --help text

## Programmatic Usage

example code

## Developer Guide

### Building and Running Correctness Tests

todo

```
1. mpicxx -std=c++17 -DCATCH_AMALGAMATED_CUSTOM_MAIN -c catch2/catch_amalgamated.cpp
2. ar rvs catch_amalgamated.a catch_amalgamated.o
```

every time:

```
mpirun -n 2 ./run_test.sh [correctness]
```

### Source Code Entry Points

todo