#ifndef __STLBFGS_H__
#define __STLBFGS_H__

#include <functional>
#include <vector>
#include <deque>
#include "linesearch.h"

#define M1QN3_PRECOND 1
namespace STLBFGS {
    typedef std::vector<double> vector;
    typedef std::function<void(const vector &x, double &f, vector &g)> func_grad_eval;
    typedef std::deque<vector> history;

    struct Optimizer {
        Optimizer(func_grad_eval func_grad) : func_grad(func_grad) {}
        bool run(vector &sol); // actual optimization loop

        struct IHessian { // L-BFGS approximates inverse Hessian matrix by storing a limited history of past updates
            void mult(const vector &g, vector &result) const; // matrix-vector multiplication
            void add_correction(const vector &s, const vector &y);

            const int history_depth;
            const bool m1qn3_precond;

            history S = {};
            history Y = {};
            vector diag = {};  // used if m1qn3_precond is set to true
            double gamma = 1.; // used otherwise
        } invH = { 10, true };

        const func_grad_eval func_grad;

        // L-BFGS user parameters
        int maxiter = 10000; // maximum number of quasi-Newton updates
        double ftol = 1e-6;  // the iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol
        double gtol = 1e-14; // the iteration stops when ||g||/max(1,||x||) <= gtol
        double gmax = 1e-14; // the iteration stops when max{|g_i|, i = 1, ..., n} <= gmax

        // Line search user parameters: the step size must satisfy Wolfe conditions with these parameters
        double mu  = 1e-4; // sufficient decrease constant (Armijo rule)
        double eta = 9e-1; // curvature condition constant
//      int lsmaxfev = 16;  // TODO move all line search parameters here

        bool verbose = true;
    };
}

#endif //__STLBFGS_H__

