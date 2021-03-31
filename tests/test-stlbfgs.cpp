#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include "stlbfgs.h"

using namespace STLBFGS;

static const double xtol = 1e-3;

TEST_CASE("2d quadratic", "[bar]") {
    std::vector<double> sol = {10, 10};

    const Optimizer::func_grad_eval func = [](std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    Optimizer opt(2, func);
    opt.run(sol);

    REQUIRE(std::abs(sol[0]-7)<xtol);
    REQUIRE(std::abs(sol[1]-1)<xtol);
}

TEST_CASE("2d Rosenbrock", "[bar]") {
    std::vector<double> sol = {-1.2, 1.0};

    const Optimizer::func_grad_eval rosenbrock = [](std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (1. - x[0])*(1. - x[0]) + 100.*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]);
        g[0] = 2.*(200.*x[0]*x[0]*x[0] - 200.*x[0]*x[1] + x[0] - 1.);
        g[1] = 200.*(x[1] - x[0]*x[0]);

    };

    Optimizer opt(2, rosenbrock);
    opt.run(sol);

    REQUIRE(std::abs(sol[0]-1)<xtol);
    REQUIRE(std::abs(sol[1]-1)<xtol);
}

