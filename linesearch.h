#ifndef __LINESEARCH_H__
#define __LINESEARCH_H__

#include <functional>

namespace STLBFGS {
    struct Sample {
        double a, f, d;
    };
    typedef std::function<Sample (const double alpha)> linesearch_function;
    void line_search(const linesearch_function phi, const Sample phi0, const double at, const double mu, const double eta);
}

#endif //__LINESEARCH_H__

