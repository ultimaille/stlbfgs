#include <catch2/catch.hpp>

#include <iostream>
#include <cmath>
#include "stlbfgs.h"

using namespace STLBFGS;
static const double xtol = 1e-3;

TEST_CASE("Powell's badly scaled function", "[L-BFGS]") {
    const Optimizer::func_grad_eval fcn = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        double f0 = pow(10, 4)*x[0]*x[1]-1;
        double f1 = exp(-x[0])+exp(-x[1])-1.0001;
        f = f0*f0+f1*f1;

        double dx0f0 = 2*pow(10, 4)*f0*x[1];
        double dx0f1 = -2*f1*exp(-x[0]);
        double dx0f = dx0f0+dx0f1;
        double dx1f0 = 2*pow(10, 4)*f0*x[0];
        double dx1f1 = -2*f1*exp(-x[1]);
        double dx1f = dx1f0+dx1f1;
        g =  {dx0f, dx1f};
    };

    std::vector<double> x = {0, 1};
    Optimizer opt(x.size(), fcn);
    opt.ftol = 0;
    opt.gtol = 0;
    opt.run(x);

    std::cerr << "Powell " << x[0] << " " << x[1] << std::endl;

    REQUIRE( std::abs(x[0]-1.098e-5) < xtol );
    REQUIRE( std::abs(x[1]-9.106) < xtol );
}

