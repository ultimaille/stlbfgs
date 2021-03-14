#ifndef __STLBFGS_H__
#define __STLBFGS_H__

#include <functional>
#include <vector>

struct LBFGSopt {
    typedef std::function<void(std::vector<double>& x, double &f, std::vector<double>& g)> func_grad_eval;

    void go(std::vector<double>& sol);

    const func_grad_eval func_grad_;
    int history_depth_ = 10;
    int max_iter_ = 10;
    bool verbose_ = false;
};

#endif //__STLBFGS_H__

