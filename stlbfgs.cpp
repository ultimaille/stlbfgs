#include <iostream>
#include "stlbfgs.h"

namespace STLBFGS {
    void LBFGSopt::optimize(std::vector<double>& x) {
        const int n = static_cast<int>(x.size());
        for (double v : x)
            std::cerr << v << std::endl;

        std::vector<double> g(n);
        double f;

        func_grad_(x, f, g);

        std::cerr << "f: " << f << std::endl << std::endl;
        for (double v : g)
            std::cerr << v << std::endl;

    }
}

