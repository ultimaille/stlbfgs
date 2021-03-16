#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include "stlbfgs.h"

typedef std::function<void(const double alpha, double &f, double &g)> func_deriv_eval;
void line_search(const func_deriv_eval phi, const double alpha0, const double mu, const double eta);
template <typename T> auto square(const T &number) { return number * number; }

TEST_CASE("foo", "[bar]") {
    // Table 1
    const func_deriv_eval func = [](const double alpha, double& f, double& g) {
        constexpr double beta = 2.;
        f = -alpha/(alpha*alpha + beta);
        g = (square(alpha)-beta)/square(beta+square(alpha));
    };
    line_search(func, 1e-3, 1e-3, 1e-1);
    REQUIRE( true );
}

