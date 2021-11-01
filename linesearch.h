#ifndef __LINESEARCH_H__
#define __LINESEARCH_H__

#include <functional>

namespace STLBFGS {
    struct Sample {
        double a, f, d; // function sample represented by the argument, the value and the derivative
    };
    typedef std::function<Sample (const double alpha)> linesearch_function;

    // xtol must be >= 0. Search fails when the relative width of the interval of uncertainty is at most xtol.
    // lsmaxfev must be > 0. Search fails when the number of calls to phi exceeds lsmaxfev.

    bool line_search(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta, const double xtol=1e-6, const int lsmaxfev=20);
}

#endif //__LINESEARCH_H__

