#include <cstdio>
#include <mpi.h>
#include <ostream>

#include <CLI11/CLI11.hpp>
#include <mtkahypar.h>

#include <types.h>
#include <hund_computation.h>
#include <kahypar_computation.h>

enum OutputType {
	PRINT,
	FILES,
};

template <typename T>
void print_vector_to_stream(std::vector<T> v, std::ostream &os) {
	for (size_t i = 0; i < v.size() - 1; i++) {
		os << v[i] << " ";
	}
	os << v[v.size() - 1] << std::endl;
}

void run_separator_size_test(
	Hypergraph hypergraph,
	BisectionConfig bisection_config,
	BreakConditionConfig break_condition_config,
	MultithreadingConfig multithreading_config
) {
	// Get the configuration right.
	BisectionConfigMtKahypar *cmk = std::get_if<BisectionConfigMtKahypar>(&bisection_config);
	BreakConditionConfigRecursionDepth *crd = std::get_if<BreakConditionConfigRecursionDepth>(&break_condition_config);
	if (!cmk || !crd) {
		printf("Error: separator-size-test needs --bisection-method to be MT_KAHYPAR and --break-condition to be RECURSION_DEPTH.\n");
		return;
	}
	int nr_of_blocks = std::pow(2, crd->depth);
	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
	if (node_id == 0) {
		printf("Running Separator Size Test\n");
		printf("Matrix dimensions: %lu x %lu, partitioned into %d blocks.\n", hypergraph.get_vertex_count(), hypergraph.get_edge_count(), nr_of_blocks);
	}

	// Run HUND algorithm.
	SeparatorSizeLogger logger = SeparatorSizeLogger();
	auto hund_computation = HUNDComputation(
		*cmk,
		*crd,
		multithreading_config,
		hypergraph,
		logger
	);
	auto hund_result = hund_computation.run_multi_node();
	logger.gather();
	auto hund_size = logger.get_total_separator_size();

	// Run KahyparComputation.
	KahyparComputation kahypar_computation(
		*cmk,
		nr_of_blocks,
		multithreading_config,
		hypergraph
	);
	auto kahypar_size = kahypar_computation.size_of_separator();

	// Print results.
	if (node_id == 0) {
		printf("HUND total size of separators: %lu\n", hund_size);
		printf("MtKaHyPar size of separator:   %lu\n", kahypar_size);
	}
}


void run_bisection_test(
	Hypergraph hypergraph,
	BisectionConfig bisection_config,
	BreakConditionConfig break_condition_config,
	MultithreadingConfig multithreading_config
) {
	// Print initial info.
	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
	if (node_id == 0) {
		printf("Running Bisection Test\n");
		printf("Matrix dimensions: %lu x %lu.\n", hypergraph.get_vertex_count(), hypergraph.get_edge_count());
	}

	// Run HUNDComputation.
	auto logger = BisectionQualityRangeLogger();
	auto hund_computation = HUNDComputation(
		bisection_config,
		break_condition_config,
		multithreading_config,
		hypergraph,
		logger
	);
	auto hund_result = hund_computation.run_multi_node();
	auto bisection_attempts = logger.gather();

	// Print recorded results.
	if (node_id == 0) {
		for (auto i = 0; i < bisection_attempts.size(); i++) {
			auto attempt = bisection_attempts[i];
			printf("Parallel Bisection #%i:\n", i);
			printf("    Recursion Depth: %i\n", attempt.recursion_depth);
			printf("    Range:           [%i, %i]\n", attempt.range_start, attempt.range_end - 1);
			printf("    Max Imbalances:  ");
			print_vector_to_stream(attempt.max_imbalances, std::cout);
			printf("    True Imbalances: ");
			print_vector_to_stream(attempt.true_imbalances, std::cout);
			printf("    Qualities:       ");
			print_vector_to_stream(attempt.qualities, std::cout);
			printf("    Rel. Sep. Sizes: ");
			print_vector_to_stream(attempt.relative_separator_sizes, std::cout);
		}
	}
}

void run_run(
	Hypergraph hypergraph,
	BisectionConfig bisection_config,
	BreakConditionConfig break_condition_config,
	MultithreadingConfig multithreading_config,
	OutputType output_type,
	std::string row_perm_file,
	std::string col_perm_file
) {
	auto logger = NoOpLogger();
	auto hund_computation = HUNDComputation(
		bisection_config,
		break_condition_config,
		multithreading_config,
		hypergraph,
		logger
	);
	auto result = hund_computation.run_multi_node();

	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
	if (node_id != 0) return;

	switch (output_type) {
	case PRINT:
		printf("Run HUND Final Result:\n");
		printf("Row permutation = ");
		print_vector(result.row_permutation);
		printf("Col permutation = ");
		print_vector(result.column_permutation);
		break;
	case FILES: {
		std::ofstream row_file;
		row_file.open(row_perm_file);
		print_vector_to_stream(result.row_permutation, row_file);
		row_file.close();

		std::ofstream col_file;
		col_file.open(col_perm_file);
		print_vector_to_stream(result.column_permutation, col_file);
		col_file.close();
		break;
	}
	default:
		assert(false);
	}
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	// 1. Configure CLI11

    CLI::App app{"HUND CLI"};
    argv = app.ensure_utf8(argv);

    int threads_per_rank;
    app.add_option("-t,--threads-per-rank", threads_per_rank, "Number of threads per MPI rank. Required. Must be at least 2.")
    	->required(true);

    std::string matrix_file = "matrix.mtx";
    app.add_option("-f,--matrix-file", matrix_file, "File path of the input matrix file.")
    	->group("Input options")
    	->required(true)
    	->check(CLI::ExistingFile);

   	std::map<std::string, MatrixFileFormat> matrix_file_format_map{{"MATRIX_MARKET", MATRIX_MARKET}};
    MatrixFileFormat matrix_file_format = MATRIX_MARKET;
    app.add_option("--matrix-file-format", matrix_file_format, "File format of the specified matrix file. Currently only supports MATRIX_MARKET.")
        ->transform(CLI::CheckedTransformer(matrix_file_format_map, CLI::ignore_case))
    	->group("Input options");

   	std::map<std::string, BisectionConfigVariant> bisection_config_variant_map{{"MT_KAHYPAR", MT_KAHYPAR}};
    BisectionConfigVariant bisection_config_variant = MT_KAHYPAR;
    app.add_option("--bisection-method", bisection_config_variant, "Method used for bisecting the hypergraph at each step. Currently only supports MT_KAHYPAR.")
        ->transform(CLI::CheckedTransformer(bisection_config_variant_map, CLI::ignore_case))
    	->group("Bisection options");

    double kahypar_max_imbalance = 0.03;
    app.add_option("--kahypar-max-imbalance", kahypar_max_imbalance, "Maximum imbalance passed as a parameter to the MT_KAHYPAR bisection method. Only has an effect if --bisection-method is set to MT_KAHYPAR.")
    	->group("Bisection options")
    	->check(CLI::PositiveNumber);

   	std::map<std::string, KahyparObjectiveFunction> kahypar_objective_function_map{{"KM1", KahyparObjectiveFunction::KM1}, {"SOED", KahyparObjectiveFunction::SOED}, {"CUT", KahyparObjectiveFunction::CUT}};
    KahyparObjectiveFunction kahypar_objective_function;
    app.add_option("--kahypar-objective-function", kahypar_objective_function, "Objective function passed as a parameter to the MT_KAHYPAR bisection method. Must be one of KM1, CUT, SOED. Only has an effect if --bisection-method is set to MT_KAHYPAR.")
        ->transform(CLI::CheckedTransformer(kahypar_objective_function_map, CLI::ignore_case))
    	->group("Bisection options");

   	std::map<std::string, BreakConditionVariant> break_condition_variant_map{{"RECURSION_DEPTH", RECURSION_DEPTH}, {"BLOCK_SIZE", BLOCK_SIZE}};
    BreakConditionVariant break_condition_variant = BLOCK_SIZE;
    app.add_option("--break-condition", break_condition_variant, "Specify at what point in the recursion to switch to a local ordering algorithm. Set to BLOCK_SIZE to stop at a certain block size (specified with --block-size), or to RECURSION_DEPTH to stop at a certain recursion depth (specified with --recursion-depth).")
        ->transform(CLI::CheckedTransformer(break_condition_variant_map, CLI::ignore_case))
    	->group("Break Condition Options");

    int block_size;
    app.add_option("--block-size", block_size, "Specify the block size at which to switch to a local ordering algorithm. A local ordering algorithm is used whenever a block's row or column count is less than or equal to --block-size. Only has an effect if --break-condition is set to BLOCK_SIZE.")
    	->group("Break Condition Options")
    	->check(CLI::NonNegativeNumber);

    int recursion_depth;
    app.add_option("--recursion-depth", recursion_depth, "Specify the recursion depth at which to switch to a local ordering algorithm. Only has an effect if --break-condition is set to RECURSION_DEPTH.")
    	->group("Break Condition Options")
    	->check(CLI::NonNegativeNumber);

    app.require_subcommand(1); // Require exactly 1 subcommand.

    CLI::App *run = app.add_subcommand("run", "Do a normal run of the HUND algorithm.");
	
   	std::map<std::string, OutputType> output_type_map{{"PRINT", PRINT}, {"FILES", FILES}};
	OutputType output_type;
    run->add_option("-o,--output-type", output_type, "Set to PRINT (default) to print the results to stdout in a human-readable way. Set to FILES to store the row and column permutation vectors in a file each, specified through --col-perm-file and --row-perm-file.")
        ->transform(CLI::CheckedTransformer(output_type_map, CLI::ignore_case))
    	->group("Output options");
    std::string col_perm_file, row_perm_file;
    run->add_option("--col-perm-file", col_perm_file, "Specify output file for the column permutation. Only has an effect if --output-type is set to FILES.")
    	->group("Output options");
    run->add_option("--row-perm-file", row_perm_file, "Specify output file for the row permutation. Only has an effect if --output-type is set to FILES.")
    	->group("Output options");

    CLI::App *bisection_test = app.add_subcommand("bisection-test", "Run the HUND algorithm as specified, and report how all the attempts at simultaneous bisection went.");

    CLI::App *separator_size_test = app.add_subcommand("separator-size-test", "Run the HUND algorithm as specified, and report the total size of separators and compare them to a run of MtKaHyPar.");

    app.footer("EXAMPLES\n"
    	"runhund \\\n"
    	"    --matrix-file=matrix.mtx --matrix-file-format=MATRIX_MARKET \\\n"
    	"    --bisection-method=MT_KAHYPAR --kahypar-max-imbalance=0.03 --kahypar-objective-function=KM1 \\\n"
    	"    --break-condition=BLOCK_SIZE --block-size=50 \\\n"
    	"    run --output-type=FILES --row-perm-file=rows.txt --col-perm-file=cols.txt"
    );

    // 2. Run CLI11.
    // Unwrap CLI11_PARSE to call MPI_Finalize(); before exiting.
    try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
    	MPI_Finalize();
    	return app.exit(e);
	}

	// 3. Create hypergraph and config objects.

	Hypergraph hypergraph = Hypergraph(matrix_file_format, matrix_file);

	MultithreadingConfig multithreading_config = {
		.number_of_threads_per_rank = threads_per_rank
	};

    BisectionConfig bisection_config;
    switch (bisection_config_variant) {
    case MT_KAHYPAR:
    	bisection_config = BisectionConfigMtKahypar {
    		.max_imbalance = kahypar_max_imbalance,
    		.objective_function = kahypar_objective_function
    	};
    	break;
    default:
    	assert(false);
    }

    BreakConditionConfig break_condition_config;
    switch (break_condition_variant) {
    case RECURSION_DEPTH:
    	break_condition_config = BreakConditionConfigRecursionDepth {
    		.depth = recursion_depth,
    	};
    	break;
    case BLOCK_SIZE:
    	break_condition_config = BreakConditionConfigBlockSize {
    		.max_block_size_inclusive = block_size,
    	};
    	break;
    default:
    	assert(false);
    }

    // 4. Run actual command.

    if (app.got_subcommand(separator_size_test)) {
    	run_separator_size_test(hypergraph, bisection_config, break_condition_config, multithreading_config);
    } else if (app.got_subcommand(bisection_test)) {
    	run_bisection_test(hypergraph, bisection_config, break_condition_config, multithreading_config);
    } else if (app.got_subcommand(run)) {
    	run_run(hypergraph, bisection_config, break_condition_config, multithreading_config, output_type, row_perm_file, col_perm_file);
    }

    MPI_Finalize();

    return 0;
}