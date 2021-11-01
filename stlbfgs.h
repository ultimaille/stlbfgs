#ifndef __STLBFGS_H__
#define __STLBFGS_H__

#include <functional>
#include <vector>
#include <deque>
#include "linesearch.h"

namespace STLBFGS {
    struct Optimizer {
        typedef std::function<void(const std::vector<double>& x, double &f, std::vector<double>& g)> func_grad_eval;
        struct IHessian { // L-BFGS approximates inverse Hessian matrix by storing a limited history of past updates
            void mult(const std::vector<double> &g, std::vector<double> &result) const; // matrix-vector multiplication
            void add_correction(const std::vector<double>&s, const std::vector<double>& y);

            int history_depth = 10;
            typedef std::deque<std::vector<double>> history;
            history S = {};
            history Y = {};
            double gamma = 1.;
        };

        void run(std::vector<double>& sol); // actual optimization loop

        const func_grad_eval func_grad;
        IHessian invH = {};  // current inverse Hessian approximation

        // L-BFGS user parameters
        int maxiter = 10000; // maximum number of quasi-Newton updates
        double ftol = 1e-6;  // the iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol
        double gtol = 1e-14; // the iteration stops when ||g||/max(1,||x||) <= gtol
        double gmax = 1e-14; // the iteration stops when max{|g_i|, i = 1, ..., n} <= gmax

        // Line search user parameters: the step size must satisfy Wolfe conditions with these parameters
        double mu  = 1e-4; // sufficient decrease constant (Armijo rule)
        double eta = 9e-1; // curvature condition constant
//      double lsxtol = 1e-6; // TODO move all line search parameters here
//      int lsmaxfev = 16;

        bool verbose = true;
    };
}

#endif //__STLBFGS_H__

