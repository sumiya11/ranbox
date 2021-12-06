#define CATCH_CONFIG_MAIN

#include "catch2/catch2.hpp"

#include "newRanBox.c"

TEST_CASE( "Simple" )
{
    int x = 1;

    CHECK( x + x == 3 );
}

TEST_CASE( "Cholesky" )
{
    vector<vector<double>> A = {
        {2, 0},
        {0, 3}
    };
    auto L = vector<vector<double>>(2, vector<double>(2, 0));
    CHECK( Cholesky( A, L ) );

    A = {
        {2, 0},
        {0, -1}
    };
    CHECK( Cholesky( A, L ) );

    A = {
        {3, 1},
        {2, 5}
    };
    CHECK( Cholesky( A, L ) );

}
