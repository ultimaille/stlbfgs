#ifndef __STLBFGS_H__
#define __STLBFGS_H__

#include <functional>
#include <vector>
#include <deque>
#include "linesearch.h"

namespace STLBFGS {
    struct Optimizer {
        typedef std::function<void(const std::vector<double>& x, double &f, std::vector<double>& g)> func_grad_eval;

        Optimizer(int nvars, func_grad_eval func) : func_grad(func), invH(nvars) {}
        void run(std::vector<double>& sol);

        const func_grad_eval func_grad;

        struct IHessian {
            IHessian(int n) : nvars(n) {}
            void mult(const std::vector<double> &g, std::vector<double> &result);
            void add_correction(const std::vector<double>&s, const std::vector<double>& y);

            const int nvars;
            int history_depth = 10;
            typedef std::deque<std::vector<double>> history;
            history S = {};
            history Y = {};
            double gamma = 1.;
        } invH;

        int maxiter = 100;   // Maximum number of function evaluations
        double ftol = 1e-6;  // The iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol
        double gtol = 1e-14; // The iteration will stop when ||g||/max(1,||x||) <= gtol
        double gmax = 1e-14; // The iteration will stop when max{|g_i|, i = 1, ..., n} <= gmax

        // Line search parameters: the step size must satisfy Wolfe conditions with these parameters
        double mu  = 1e-4; // sufficient decrease constant (Armijo rule)
        double eta = 9e-1; // curvature condition constant, TODO try 1e-2
        double stpmin = 1e-15;
        double stpmax = 1e15;
        int lsmaxfev = 20;

        bool verbose = true;
    };
}

#endif //__STLBFGS_H__

