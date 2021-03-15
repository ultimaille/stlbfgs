#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include "stlbfgs.h"

static const double xtol = 1e-3;

TEST_CASE("foo", "[bar]") {
    std::vector<double> sol = {10, 10};

    const LBFGSopt::func_grad_eval func = [](std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    LBFGSopt opt(2, func);
    opt.go(sol);

    REQUIRE(std::abs(sol[0]-7)<xtol);
    REQUIRE(std::abs(sol[1]-1)<xtol);
}

