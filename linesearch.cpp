#include <iostream>
#include <functional>
//#include <vector>
#include <cmath>
#include <limits>
#include <cassert>

/*
typedef std::vector<double> vector;
typedef std::function<void(std::vector<double>& x, double &f, std::vector<double>& g)> func_grad_eval;


void line_search(func_grad_eval phi,
        double step0,
        double alpha = 1.,
        double c1 = 1e-4,
        double c2 = 1e-1,
        double xtol = 1e-14,
        double alpha_min = 0.,
        double alpha_max = std::numeric_limits<double>::max(),
        int maxfev = std::numeric_limits<int>::max(),
        double delta = 0.66) {

        constexpr double xtrapf = 4.;
        if (alpha <= 0) {
            std::cerr << "alpha < 0.0" << std::endl;
            return;
        }
        if (c1 < 0) {
            std::cerr << "c1 < 0.0" << std::endl;
            return;
        }
        if (c2 < 0) {
            std::cerr << "c2 < 0.0" << std::endl;
            return;
        }
        if (xtol < 0) {
            std::cerr << "xtol < 0.0" << std::endl;
            return;
        }
        if (alpha_min < 0) {
            std::cerr << "alpha_min < 0.0" << std::endl;
            return;
        }
        if (alpha_max < 0) {
            std::cerr << "alpha_max < 0.0" << std::endl;
            return;
        }
        if (maxfev <= 0) {
            std::cerr << "maxfev < 0" << std::endl;
            return;
        }
}
*/

