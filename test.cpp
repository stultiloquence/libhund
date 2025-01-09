#include <mpi.h>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_amalgamated.hpp>

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

SCENARIO("Sequential Testing", "[2rank]") {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}




int main( int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    return result;
}