#define CATCH_CONFIG_RUNNER
#include <catch2/catch_amalgamated.hpp>
#include <mpi.h>

#include <hund_computation.h>
#include <kahypar_computation.h>


void require_permutation(std::vector<unsigned long> vector) {
    std::sort(vector.begin(), vector.end());
    for (int i = 0; i < vector.size(); i++) {
        REQUIRE(i == vector[i]);
    }
}

TEST_CASE("Obvious block structure is identified", "[correctness]") {
    std::vector<size_t> example_hyperedge_indices { 0, 2, 4, 6, 8, 10, 12, 14, 16 };
    std::vector<unsigned long> example_hyperedges { 0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7 };

    auto hypergraph = Hypergraph(8, example_hyperedge_indices, example_hyperedges);
    auto bisector = MtKahyparBisector(
        0.05,
        0
    );
    auto computation = HUNDComputation(
        bisector,
        (BreakConditionConfigBlockSize) {
            .max_block_size_inclusive = 2,
        },
        hypergraph
    );

    auto result = computation.run_multi_node();

    auto row_permutation = result.row_permutation;
    auto column_permutation = result.column_permutation;
    require_permutation(row_permutation);
    require_permutation(column_permutation);
    for (int i = 0; i < 8; i += 2) {
        long ri = (long) row_permutation[i], ri1 = (long) row_permutation[i + 1];
        long ci = (long) column_permutation[i], ci1 = (long) column_permutation[i + 1];
        REQUIRE(std::abs(ri - ri1) == 1);
        REQUIRE(std::abs(ci - ci1) == 1);
    }
}

TEST_CASE("Time to run on real life matrix", "[performance]") {
    const int recursion_depth = 3;
    const int nr_of_blocks = pow(2, recursion_depth);
    const double max_imbalance = 0.05;

    auto hypergraph = Hypergraph(MATRIX_MARKET, "../examples/fd18.mtx");

    auto bisector = MtKahyparBisector(max_imbalance);
    auto computation = HUNDComputation(
        bisector,
        (BreakConditionConfigRecursionDepth) {
            .depth = recursion_depth,
        },
        hypergraph
    );

    BENCHMARK("fd18 HUNDComputation.run_multi_node()") {
        return computation.run_multi_node();
    };

    KahyparComputation kahypar_computation(
        max_imbalance,
        0,
        nr_of_blocks,
        hypergraph
    );

    BENCHMARK("fd18 KahyparComputation.dissect()") {
        return kahypar_computation.dissect();
    };
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    return result;
}