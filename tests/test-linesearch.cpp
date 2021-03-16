#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>

#include "stlbfgs.h"

double find_cubic_minimizer(double a, double fa, double ga, double b, double fb, double gb);
double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb);
double find_quadratic_minimizer(double a, double ga, double b, double gb);

TEST_CASE("cubic 1", "[interpolation]") {
    auto f = [](const double a) { return a*a*a - a; };
    auto g = [](const double a) { return 3*a*a - 1; };
    double a = -.4;
    double b = 19.;
    double minimizer = find_cubic_minimizer(a, f(a), g(a), b, f(b), g(b));
    REQUIRE( std::abs(minimizer-1./std::sqrt(3.))<1e-14 );
}

TEST_CASE("cubic 2", "[interpolation]") {
    auto f = [](const double a) { return a*a*a - 3*a*a - 144*a + 432; };
    auto g = [](const double a) { return 3*a*a - 6*a - 144; };
    double a = -5.;
    double b = 113.;
    double minimizer = find_cubic_minimizer(a, f(a), g(a), b, f(b), g(b));
    REQUIRE( std::abs(minimizer-8.)<1e-14 );
}

TEST_CASE("quadratic", "[interpolation]") {
    for (int i=0; i<8; i++) {
        double A = rand()/(double)RAND_MAX + .2134;
        double B = rand()/(double)RAND_MAX;
        double C = rand()/(double)RAND_MAX;
        double M = -B/(2.*A);
        auto f = [&](const double a) { return A*a*a + B*a + C; };
        auto g = [&](const double a) { return 2*A*a + B; };
        double a = M - rand()/(double)RAND_MAX;
        double b = M + rand()/(double)RAND_MAX;
        double minimizer1 = find_quadratic_minimizer(a, f(a), g(a), b, f(b));
        REQUIRE( std::abs(minimizer1 - M)<1e-14 );
        double minimizer2 = find_quadratic_minimizer(a, g(a), b, g(b));
        REQUIRE( std::abs(minimizer2 - M)<1e-14 );
    }
}

template <typename T> auto square(const T &number) { return number * number; }
typedef std::function<void(const double alpha, double &f, double &g)> func_deriv_eval;
void line_search(const func_deriv_eval phi, const double alpha0, const double mu, const double eta);


/*
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
*/

