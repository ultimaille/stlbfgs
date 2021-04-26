#include <iostream>
#include <deque>
#define _USE_MATH_DEFINES
#include <cmath>
#undef NDEBUG
#include <cassert>
#include "stlbfgs.h"

namespace STLBFGS {

    typedef std::vector<double> vector;

    // compute dot product <a,b>
    double dot(const vector &a, const vector &b) {
        assert(a.size()==b.size());
        double dot = 0;
#pragma omp parallel for reduction(+:dot)
        for (int i=0; i<(int)a.size(); i++)
            dot += a[i]*b[i];
        return dot;
    }

    double norm(const vector &v) {
        return std::sqrt(dot(v, v));
    }

    // Add a correction pair {s, y} to the optimization history
    void Optimizer::IHessian::add_correction(const vector &s, const vector &y) {
        assert(nvars == (int)s.size());
        assert(nvars == (int)y.size());
        S.push_front(s);
        Y.push_front(y);
        if ((int)S.size()>history_depth) S.pop_back();
        if ((int)Y.size()>history_depth) Y.pop_back();

        double yy = dot(y, y);
        assert(std::abs(yy) > 0);
        gamma = dot(s, y)/yy;
        assert(std::isfinite(gamma));
    }

    // Multiply a vector g by the inverse Hessian matrix approximation
    // Algorithm 7.4 (L-BFGS two-loop recursion)
    // Nocedal and Wright, Numerical optimization (2006)
    void Optimizer::IHessian::mult(const vector &g, vector &result) {
        const int m = static_cast<int>(S.size());
        assert((int)Y.size() == m);
        assert((int)g.size() == nvars);

        result = g;

        std::vector<double> a(m);
        for (int i=0; i<m; i++) {
            const vector &y = Y[i];
            const vector &s = S[i];
            double sy = dot(s, y);
            assert(std::abs(sy) > 0);
            a[i] = dot(s, result)/sy;
            assert(std::isfinite(a[i]));
#pragma omp parallel for
            for (int j=0; j<nvars; j++)
                result[j] -= a[i]*y[j];
        }

#pragma omp parallel for
        for (int j=0; j<nvars; j++)
            result[j] *= gamma;

        for (int i=m; i--;) {
            const vector &y = Y[i];
            const vector &s = S[i];
            double b = dot(y, result)/dot(s, y);
            assert(std::isfinite(b));
#pragma omp parallel for
            for (int j=0; j<nvars; j++)
                result[j] += (a[i]-b)*s[j];
        }
    }

    // Run an L-BFGS optimization
    // Algorithm 7.5
    // Nocedal and Wright, Numerical optimization (2006)
    void Optimizer::run(vector &x) {
        const int n = (int)x.size();
        assert(invH.nvars == n);

        double f;
        vector g(n), p(n);

        func_grad(x, f, g);
        for (int i=0; i<maxiter; i++) {
            invH.mult(g, p);

            assert(-dot(g, p)<0);

            double fprev = f;
            vector xprev = x, gprev = g;

            int nfev = 0;
            const linesearch_function ls_func = [&](const double alpha) -> Sample {
                nfev++;
#pragma omp parallel for
                for (int j=0; j<n; j++)
                    x[j] = xprev[j]-p[j]*alpha;
                func_grad(x, f, g);
                return { alpha, f, -dot(g, p) };
            };

            double alpha = i ? 1. : 1./norm(g);
            assert(std::isfinite(alpha));
            bool res = line_search(ls_func, {0, f, -dot(g, p)}, alpha, mu, eta);
            if (!res) {
                if (verbose) std::cerr << "Linesearch failed, step = " << alpha << std::endl;
                break;
            }

            vector s(n), y(n);
#pragma omp parallel for
            for (int j=0; j<n; j++) {
                s[j] = x[j]-xprev[j];
                y[j] = g[j]-gprev[j];
            }
            invH.add_correction(s, y);

            if ((fprev-f)/std::max(std::max(std::abs(fprev), std::abs(f)), 1.) <= ftol) {
                if (verbose) std::cerr << "(f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= " << ftol << std::endl;
                break;
            }

            if (norm(g)/std::max(1., norm(x))<=gtol) {
                if (verbose) std::cerr << "||g||/max(1,||x||) <= " << gtol << std::endl;
                break;
            }

            double gmax_ = 0.;
#if defined(_OPENMP) && _OPENMP>=200805
#pragma omp parallel for reduction(max:gmax_)
#endif
            for (double gi : g)
                gmax_ = std::max(gmax_, std::abs(gi));
            if (gmax_ <= gmax) {
                if (verbose) std::cerr << "max{|g_i|, i = 1, ..., n} <= " <<  gmax << std::endl;
                break;
            }

            if (i==maxiter-1) {
                if (verbose) std::cerr << "reached maxiter" << std::endl;
            }
        }
    }

}

