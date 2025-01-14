#include <cstdio>
#include <mpi.h>
#include <ostream>

#include <include/CLI11.hpp>
#include <mtkahypar.h>

#include <types.h>
#include <hund_computation.h>
#include <kahypar_computation.h>

template <typename T>
void print_vector_to_stream(std::vector<T> v, std::ostream &os) {
	for (size_t i = 0; i < v.size() - 1; i++) {
		os << v[i] << " ";
	}
	os << v[v.size() - 1] << std::endl;
}

void run_separator_size_test(std::string &matrix_file) {
	const int recursion_depth = 3;
	const int nr_of_blocks = pow(2, recursion_depth);

	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
	if (node_id == 0) {
		printf("Running Separator Size Test\n");
	}

	SeparatorSizeLogger logger = SeparatorSizeLogger();
	auto hypergraph = Hypergraph(MATRIX_MARKET, matrix_file);
	BisectionConfigMtKahypar kahypar_config {
		.max_imbalance = 0.03,
		.objective_function = HUND_KM1,
	};

	if (node_id == 0) {
		printf("Matrix dimensions: %lu x %lu, partitioned in %d blocks.\n", hypergraph.get_vertex_count(), hypergraph.get_edge_count(), nr_of_blocks);
	}

	auto hund_computation = HUNDComputation(
		kahypar_config,
		(BreakConditionConfigRecursionDepth) {
			.depth = recursion_depth,
		},
		hypergraph,
		logger
	);

	auto hund_result = hund_computation.run_multi_node();
	logger.gather();
	auto hund_size = logger.get_total_separator_size();

	KahyparComputation kahypar_computation(kahypar_config, nr_of_blocks, hypergraph);
	auto kahypar_size = kahypar_computation.size_of_separator();

	if (node_id == 0) {
		printf("HUND total size of separators: %lu\n", hund_size);
		printf("MtKaHyPar size of separator:   %lu\n", kahypar_size);
	}
}


void run_bisection_test(std::string &matrix_file) {
	const int recursion_depth = 3;
	const int nr_of_blocks = pow(2, recursion_depth);

	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
	if (node_id == 0) {
		printf("Running Bisection Test\n");
	}

	auto logger = BisectionQualityRangeLogger();
	auto hypergraph = Hypergraph(MATRIX_MARKET, matrix_file);
	BisectionConfigMtKahypar kahypar_config {
		.max_imbalance = 0.03,
		.objective_function = HUND_KM1,
	};

	if (node_id == 0) {
		printf("Matrix dimensions: %lu x %lu, partitioned in %d blocks.\n", hypergraph.get_vertex_count(), hypergraph.get_edge_count(), nr_of_blocks);
	}

	auto hund_computation = HUNDComputation(
		kahypar_config,
		(BreakConditionConfigRecursionDepth) {
			.depth = recursion_depth,
		},
		hypergraph,
		logger
	);

	auto hund_result = hund_computation.run_multi_node();
	auto bisection_attempts = logger.gather();

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

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

    CLI::App app{"HUND CLI"};
    argv = app.ensure_utf8(argv);

    std::string matrix_file = "matrix.mtx";
    app.add_option("-f,--matrix-file", matrix_file, "File path of the input matrix file.");

    CLI::App *bisection_test = app.add_subcommand("bisection-test", "Report how all the attempts at simultaneous bisection went.");

    CLI::App *separator_size_test = app.add_subcommand("separator-size-test", "Report the total size of separators and compare them to a run of MtKaHyPar.");

    CLI11_PARSE(app, argc, argv);

    if (app.got_subcommand(separator_size_test)) {
    	run_separator_size_test(matrix_file);
    }

    if (app.got_subcommand(bisection_test)) {
    	run_bisection_test(matrix_file);
    }

    MPI_Finalize();

    return 0;
}

void old_main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	auto hypergraph = Hypergraph(MATRIX_MARKET, argv[1]);

	SeparatorSizeLogger logger = SeparatorSizeLogger();
	auto computation = HUNDComputation(
		(BisectionConfigMtKahypar) {
			.max_imbalance = 0.03,
			.objective_function = HUND_KM1,
		},
		(BreakConditionConfigRecursionDepth) {
			.depth = 3,
		},
		hypergraph,
		logger
	);

	auto result = computation.run_multi_node();

	int node_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);

	if (node_id == 0) {
		// printf("Run HUND Final Result:\n");
		// printf("Row permutation = ");
		// print_vector(result.row_permutation);
		// printf("Col permutation = ");
		// print_vector(result.column_permutation);


		std::ofstream row_file;
		row_file.open(argv[2]);
		print_vector_to_stream(result.row_permutation, row_file);
		row_file.close();

		std::ofstream col_file;
		col_file.open(argv[3]);
		print_vector_to_stream(result.column_permutation, col_file);
		col_file.close();
	}

	MPI_Finalize();
}
