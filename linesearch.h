#ifndef __LINESEARCH_H__
#define __LINESEARCH_H__

#include <functional>

namespace STLBFGS {
    struct Sample {
        double a, f, d;
    };
    typedef std::function<Sample (const double alpha)> linesearch_function;
    void line_search(const linesearch_function phi, double at, double mu, double eta);
}

#endif //__LINESEARCH_H__

