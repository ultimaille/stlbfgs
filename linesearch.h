#ifndef __LINESEARCH_H__
#define __LINESEARCH_H__

#include <functional>

namespace STLBFGS {
    struct Sample {
        double a, f, d;
    };
    typedef std::function<Sample (const double alpha)> linesearch_function;
    bool line_search(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta, const double stpmin=1e-15, const double stpmax=1e15, const int lsmaxfev=20);
}

#endif //__LINESEARCH_H__

