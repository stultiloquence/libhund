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

Example Usage:

```bash
hundcli \
    --matrix-file=examples/paper_example.mtx --matrix-file-format=MATRIX_MARKET \
    --bisection-method=MT_KAHYPAR --kahypar-max-imbalance=0.05 \
    --break-condition=BLOCK_SIZE --block-size=2 \
    run --output-type=FILES --row-file=results/rows.txt --col-file=result/cols.txt
```

Detailed Options:

```
HUND CLI

Usage: hundcli [General Options] [Input Options] [Bisection Options]
		[Break Condition Options] SUBCOMMAND [Subcommand Options]



General Options:

-h,--help
		Print this help message and exit



Input Options:
 
-f,--matrix-file TEXT:FILE REQUIRED
		File path of the input matrix file.
 
--matrix-file-format ENUM:value in {MATRIX_MARKET->0} OR {0}
		File format of the specified matrix file. Currently only supports
		MATRIX_MARKET.



Bisection Options:
 
--bisection-method ENUM:value in {MT_KAHYPAR->0} OR {0}
		Method used for bisecting the hypergraph at each step. Currently only
		supports MT_KAHYPAR.
 
--kahypar-max-imbalance FLOAT:POSITIVE
		Maximum imbalance passed as a parameter to the MT_KAHYPAR bisection
		method. Only has an effect if --bisection-method is set to MT_KAHYPAR.
		Default value is 0.03.
 
--kahypar-threads INT:NONNEGATIVE
		Number of threads per MPI rank when using MtKaHyPar. Default is 0, which
		means use the maximum number available.



Break Condition Options:
 
--break-condition ENUM:value in {BLOCK_SIZE->1,RECURSION_DEPTH->0} OR {1,0}
		Specify at what point in the recursion to switch to a local ordering
		algorithm. Set to BLOCK_SIZE to stop at a certain block size (specified
		with --block-size), or to RECURSION_DEPTH to stop at a certain recursion
		depth (specified with --recursion-depth). Default value is BLOCK_SIZE.
 
--block-size INT:NONNEGATIVE
		Specify the block size at which to switch to a local ordering algorithm.
		A local ordering algorithm is used whenever a block's row or column
		count is less than or equal to --block-size. Only has an effect if
		--break-condition is set to BLOCK_SIZE. Default value is 10.
 
--recursion-depth INT:NONNEGATIVE
		Specify the recursion depth at which to switch to a local ordering
		algorithm. Only has an effect if --break-condition is set to
		RECURSION_DEPTH.  Default value is 1.



Subcommands and their Options:

run
		Do a normal run of the HUND algorithm.

		-h,--help
		        Print this help message and exit

		-o,--output-type ENUM:value in {FILES->1,PRINT->0} OR {1,0}
				Set to PRINT to print the results to stdout in a human-readable
				way. Set to FILES to store the row and column permutation
				vectors in a file each, specified through --col-file and
				--row-file. Default value is PRINT.

         --col-file TEXT
         		Specify output file for the column permutation. Only has an
         		effect if --output-type is set to FILES. Default value is
         		col_perms.txt.

		--row-file TEXT
				Specify output file for the row permutation. Only has an effect
				if --output-type is set to FILES. Default value is
				row_perms.txt.

bisection-test

		Run the HUND algorithm as specified, and report how all the attempts at
		simultaneous bisection went.

		-h,--help
				Print this help message and exit

		--output-file TEXT
				Output file for the bisetion test result (in JSON format).
				Default value is bisection-test_results.json

separator-size-test

		Run the HUND algorithm as specified, and report the total size of
		separators and compare them to a run of MtKaHyPar.

		-h,--help
				Print this help message and exit
```

## Programmatic Usage

The three examples below should give a good idea on how to use libhund programmatically. For a full description of all methods and parameters, see the documentation in the source code.

### Example 1: Configuring Hypergraph Manually

```CPP
size_t vertex_count = 8;
std::vector<size_t> example_hyperedge_indices { 0, 2, 4, 6, 8, 10, 12, 14, 16 };
std::vector<unsigned long> example_hyperedges { 0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7 };

auto hypergraph = Hypergraph(vertex_count, example_hyperedge_indices, example_hyperedges);
auto bisector = MtKahyparBisector(
    0.05, // maximal imbalance
    0 // number of threads. 0 = use all available.
);
auto break_condition = (BreakConditionConfigBlockSize) {
	.max_block_size_inclusive = 2,
};

auto computation = HUNDComputation(
    bisector,
    break_condition,
    hypergraph
);
auto result = computation.run_multi_node();

std::count << "Row Permutation:" << std::endl;
for (auto v : result.row_permutation) {
	std::cout << v << ' ';
}
std::count << std::endl << "Column Permutation:" << std::endl;
for (auto v : result.column_permutation) {
	std::cout << v << ' ';
}
```

### Example 2: Loading Hypergraph from File

```CPP
auto hypergraph = Hypergraph(
	MATRIX_MARKET, // Only supported file format at the time.
	"examples/fd18.mtx"
);
auto bisector = MtKahyparBisector(0.05, 0);
auto break_condition = (BreakConditionConfigRecursionDepth) {
	.depth = 4,
};

auto computation = HUNDComputation(
    bisector,
    break_condition,
    hypergraph
);
auto result = computation.run_multi_node();
//...
```

### Example 3: Using a Logger

```CPP
auto hypergraph = Hypergraph(MATRIX_MARKET, "examples/fd18.mtx");
auto bisector = MtKahyparBisector(0.05, 0);
auto break_condition = (BreakConditionConfigRecursionDepth) { .depth = 4 };

auto logger = BisectionQualityRangeLogger(); // Or any other Logger from logger.h. It is also possible to write your own logger by deriving from the Logger bsae class.
auto computation = HUNDComputation(
    bisector,
    break_condition,
    hypergraph,
    logger
);
auto result = computation.run_multi_node();

// Gathering and retrieving the logger results depends on the specific logger.
auto logging_result = logger.gather();
```

## Developer Guide

### Building and Running Correctness Tests

In order to build the correctness tests, the Catch2 library file has to be built once as follows:

```bash
cd external/catch2
mpicxx -std=c++17 -DCATCH_AMALGAMATED_CUSTOM_MAIN -c catch_amalgamated.cpp
ar rvs catch_amalgamated.a catch_amalgamated.o
rm catch_amalgamated.o
```

To build the tests themselves, you then simply run (from the libhund root folder):

```bash
test/build_test.sh
```

This creates the `test/test` binary. To run the tests, you can either run the test binary directly (using a single mpi node), like so:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(realpath ./external/mt-kahypar/build/lib/) 
cd test
./test [correctness]
```

The `[correctness]` argument specifies the type of tests to be run. The two available options are `[correctness]` and `[performance]`. Running the performance tests can take a couple of minutes.

Or you can run the tests on multiple MPI nodes with `mpiexec`. To suppress output from all but one node, the `test/run_test.sh` wrapper script is provided (which also automatically changes into the `test` directory):

```bash
mpiexec -n 4 ./test/run_test.sh [correctness]
```

### Source Code Entry Points

HUND is a header-only library. Here is a brief overview over the individual files' responsibilities:

- `hypergraph.h`: Provides the `Hypergraph` class, which represents hypergraphs in compressed column storage. Hypergraphs can either be loaded from file or constructed manually through passing vectors.
- `bisector.h`: Provides `Bisector`s, which represent a specific hypergraph bisection method to be used in a single step of the HUND algorithm. The only algorithm currently available is the `MtKahyparBisector`. The derived classes need to implement the `bisect` method, which only needs to return the vertex partitioning of the bisection.
- `hypergraph_bisection.h`: The `HypergraphBisection` class is constructed from a `Hypergraph` instance together with vertex partition (usually obtained from a `Bisector`), and computes the associated sub-hypergraphs, edge partition and vertex and edge permutations.
- `hund_computation.h`: An instance of the `HUNDComputation` class represents a single fully configured run of the HUND algorithm. The configuration happens through its constructor(s), and the run can be executed through the `run_single_node` or `run_multi_node` methods. Internally, the class utilizes the `Bisector` and `HypergraphBisection` classes.
- `logger.h`: Provides the abstract `Logger` class together with some concrete loggers such as `NoOpLogger`, `SeparatorSizeLogger` or `BisectionQualityRangeLogger`. Logggers can be passed to a `HUNDComputation` and provide hooks that allow recording various information about the computation.
- `kahypar_computation.h`: The `KahyparComputation` class is only used for the `separator-size-test` subcommand, which compares the separator sizes of the `HUNDComputation` with that of the `KahyparComputation`. Both use the `MtKaHyPar` library for hypergraph dissection, but while `HUNDComputation` bisects recursively, `KahyparComputation` directly asks the library to dissect the hypergraph into the desired final number of blocks.
- `initialize_mt_kahypar.h`: Helper file providing the `initialize_mt_kahypar()` function, that prevents initializing the MtKaHyPar library more than once.
- `util.h`: Provides functions to print vectors, and functions to compute simple operations on permutations that are used by the `HUNDComputation`.
- `hund_cli.cpp`: Puts all of the above together to provide the CLI. Uses the `CLI11` library.
- `test/test.cpp`: Implements the `Catch2` tests (with a custom main function). Has correctness and performance tests.