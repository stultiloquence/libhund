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

enum BisectionConfigVariant { 
	MT_KAHYPAR 
};


enum BreakConditionVariant {
  RECURSION_DEPTH,
  BLOCK_SIZE,
};

void run_separator_size_test(
	Hypergraph hypergraph,
	Bisector &bisector,
	BreakConditionConfig break_condition_config,
	BisectionConfigVariant bisection_config_variant,
	double kahypar_max_imbalance,
	int threads_per_rank
) {
	// Get the configuration right.
	BreakConditionConfigRecursionDepth *crd = std::get_if<BreakConditionConfigRecursionDepth>(&break_condition_config);
	if (bisection_config_variant != MT_KAHYPAR || !crd) {
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
	SeparatorSizeLogger logger = SeparatorSizeLogger(
		crd->depth
	);
	auto hund_computation = HUNDComputation(
		bisector,
		*crd,
		hypergraph,
		logger
	);
	auto hund_result = hund_computation.run_multi_node();
	logger.gather();
	auto hund_size = logger.get_total_separator_size();
	auto hund_size_weighted = logger.get_total_separator_size_weighted();

	// Run KahyparComputation.
	KahyparComputation kahypar_computation(
		kahypar_max_imbalance,
		threads_per_rank,
		nr_of_blocks,
		hypergraph
	);
	auto kahypar_size = kahypar_computation.size_of_separator();
	auto kahypar_cut = kahypar_computation.get_cut();
	auto kahypar_km1 = kahypar_computation.get_km1();
	auto kahypar_soed = kahypar_computation.get_soed();

	// Print results.
	if (node_id == 0) {
		printf("HUND total size of separators:          %lu\n", hund_size);
		printf("HUND total size of separators weighted: %lu\n", hund_size_weighted);
		printf("MtKaHyPar size of separator:            %lu\n", kahypar_size);
		printf("MtKaHyPar CUT:                          %f\n", kahypar_cut);
		printf("MtKaHyPar KM1:                          %f\n", kahypar_km1);
		printf("MtKaHyPar SOED:                         %f\n", kahypar_soed);
	}
}


void run_bisection_test(
	Hypergraph hypergraph,
	Bisector &bisector,
	BreakConditionConfig break_condition_config,
	std::string output_file_name
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
		bisector,
		break_condition_config,
		hypergraph,
		logger
	);
	auto hund_result = hund_computation.run_multi_node();
	auto bisection_attempts = logger.gather();

	// Print recorded results.
	if (node_id == 0) {
		std::ofstream output_file;
		output_file.open(output_file_name);
		print_as_json(bisection_attempts, output_file);
		output_file.close();
	}
}

void run_run(
	Hypergraph hypergraph,
	Bisector &bisector,
	BreakConditionConfig break_condition_config,
	OutputType output_type,
	std::string row_perm_file,
	std::string col_perm_file
) {
	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);

	if (node_id == 0) {
		printf("Run HUND Computation\n");
		printf("Matrix dimensions: %lu x %lu\n", hypergraph.get_vertex_count(), hypergraph.get_edge_count());
	}
	auto hund_computation = HUNDComputation(
		bisector,
		break_condition_config,
		hypergraph
	);
	auto result = hund_computation.run_multi_node();

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

    int threads_per_rank = 0;
    app.add_option("-t,--threads-per-rank", threads_per_rank, "Number of threads per MPI rank. Default is 0, which means use the maximum number available.")
    	->check(CLI::NonNegativeNumber);

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
	app.add_option("--kahypar-max-imbalance", kahypar_max_imbalance, "Maximum imbalance passed as a parameter to the MT_KAHYPAR bisection method. Only has an effect if --bisection-method is set to MT_KAHYPAR. Default value is 0.03.")
    	->group("Bisection options")
    	->check(CLI::PositiveNumber);

   	std::map<std::string, BreakConditionVariant> break_condition_variant_map{{"RECURSION_DEPTH", RECURSION_DEPTH}, {"BLOCK_SIZE", BLOCK_SIZE}};
    BreakConditionVariant break_condition_variant = BLOCK_SIZE;
    app.add_option("--break-condition", break_condition_variant, "Specify at what point in the recursion to switch to a local ordering algorithm. Set to BLOCK_SIZE to stop at a certain block size (specified with --block-size), or to RECURSION_DEPTH to stop at a certain recursion depth (specified with --recursion-depth). Default value is BLOCK_SIZE.")
        ->transform(CLI::CheckedTransformer(break_condition_variant_map, CLI::ignore_case))
    	->group("Break Condition Options");

    int block_size = 10;
    app.add_option("--block-size", block_size, "Specify the block size at which to switch to a local ordering algorithm. A local ordering algorithm is used whenever a block's row or column count is less than or equal to --block-size. Only has an effect if --break-condition is set to BLOCK_SIZE. Default value is 10.")
    	->group("Break Condition Options")
    	->check(CLI::NonNegativeNumber);

    int recursion_depth = 1;
    app.add_option("--recursion-depth", recursion_depth, "Specify the recursion depth at which to switch to a local ordering algorithm. Only has an effect if --break-condition is set to RECURSION_DEPTH.  Default value is 1.")
    	->group("Break Condition Options")
    	->check(CLI::NonNegativeNumber);

    app.require_subcommand(1); // Require exactly 1 subcommand.

    CLI::App *run = app.add_subcommand("run", "Do a normal run of the HUND algorithm.");
	
   	std::map<std::string, OutputType> output_type_map{{"PRINT", PRINT}, {"FILES", FILES}};
	OutputType output_type = PRINT;
    run->add_option("-o,--output-type", output_type, "Set to PRINT to print the results to stdout in a human-readable way. Set to FILES to store the row and column permutation vectors in a file each, specified through --col-perm-file and --row-perm-file. . Default value is PRINT.")
        ->transform(CLI::CheckedTransformer(output_type_map, CLI::ignore_case))
    	->group("Output options");
    std::string col_perm_file = "col_perms.txt";
	std::string row_perm_file = "row_perms.txt";
    run->add_option("--col-perm-file", col_perm_file, "Specify output file for the column permutation. Only has an effect if --output-type is set to FILES. Default value is col_perms.txt.")
    	->group("Output options");
    run->add_option("--row-perm-file", row_perm_file, "Specify output file for the row permutation. Only has an effect if --output-type is set to FILES. Default value is row_perms.txt.")
    	->group("Output options");

    CLI::App *bisection_test = app.add_subcommand("bisection-test", "Run the HUND algorithm as specified, and report how all the attempts at simultaneous bisection went.");

    std::string output_file = "bisection-test_results.json";
    bisection_test->add_option("--output-file", output_file, "Output file for the bisetion test result (in JSON format). Default value is bisection-test_results.json");

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

    std::unique_ptr<Bisector> bisector;
    switch (bisection_config_variant) {
    case MT_KAHYPAR:
    	bisector = std::make_unique<MtKahyparBisector>(
    		kahypar_max_imbalance, 
    		threads_per_rank
    	);
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
    	run_separator_size_test(hypergraph, *bisector, break_condition_config, bisection_config_variant, kahypar_max_imbalance, threads_per_rank);
    } else if (app.got_subcommand(bisection_test)) {
    	run_bisection_test(hypergraph, *bisector, break_condition_config, output_file);
    } else if (app.got_subcommand(run)) {
    	run_run(hypergraph, *bisector, break_condition_config, output_type, row_perm_file, col_perm_file);
    }

    MPI_Finalize();

    return 0;
}