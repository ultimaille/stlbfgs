#include <iostream>
#include <deque>
#include <cmath>
#include <cassert>
#include "stlbfgs.h"

namespace STLBFGS {

    typedef std::vector<double> vector;

    // compute dot product <a,b>
    double operator*(const vector &a, const vector &b) {
        assert(a.size()==b.size());
        double dot = 0;
#pragma omp parallel for reduction(+:dot)
        for (size_t i=0; i<a.size(); i++)
            dot += a[i]*b[i];
        return dot;
    }

    double norm(const vector &v) {
        return std::sqrt(v*v);
    }

    // Add a correction pair {s, y} to the optimization history
    void Optimizer::IHessian::add_correction(const vector &s, const vector &y) {
        assert(nvars == s.size());
        assert(nvars == y.size());
        S.push_front(s);
        Y.push_front(y);
        if (S.size()>history_depth) S.pop_back();
        if (Y.size()>history_depth) Y.pop_back();

        double yy = y*y;
        assert(std::abs(yy) > 1e-14);
        gamma = (s*y)/yy;
    }

    // Multiply a vector g by the inverse Hessian matrix approximation
    // Algorithm 7.4 (L-BFGS two-loop recursion)
    // Nocedal and Wright, Numerical optimization (2006)
    void Optimizer::IHessian::mult(const vector &g, vector &result) {
        const size_t m = S.size();
        assert(Y.size() == m);
        assert(g.size() == nvars);

        result = g;

        std::vector<double> a(m);
        for (size_t i=0; i<m; i++) {
            const vector &y = Y[i];
            const vector &s = S[i];
            double sy = s*y;
            assert(std::abs(sy) > 1e-14);
            a[i] = (s*result)/sy;
#pragma omp parallel for
            for (size_t j=0; j<nvars; j++)
                result[j] -= a[i]*y[j];
        }

        if (m>0) {
#pragma omp parallel for
            for (size_t j=0; j<nvars; j++)
                result[j] *= gamma;
        }

        for (size_t i=m; i--;) {
            const vector &y = Y[i];
            const vector &s = S[i];
            double b = (y*result)/(s*y);
#pragma omp parallel for
            for (size_t j=0; j<nvars; j++)
                result[j] += (a[i]-b)*s[j];
        }
    }

    // Run an L-BFGS optimization
    // Algorithm 7.5
    // Nocedal and Wright, Numerical optimization (2006)
    void Optimizer::run(vector &x) {
        const size_t n = x.size();
        assert(invH.nvars == n);

        double f;
        vector g(n), p(n);

        func_grad(x, f, g);
        for (size_t i=0; i<maxiter; i++) {
            {
                std::cerr << "x: ";
                for (double v : x) std::cerr << v << " ";
                std::cerr << std::endl;
                std::cerr << "g: ";
                for (double v : g) std::cerr << v << " ";
                std::cerr << std::endl;
                std::cerr << "f: " << f << std::endl;
            }

            invH.mult(g, p);

            double fprev = f;
            vector xprev = x;
            vector gprev = g;
            vector s(n), y(n);

            for (size_t j=0; j<n; j++)
                x[j] -= p[j]/3.;
            func_grad(x, f, g);

#pragma omp parallel for
            for (size_t j=0; j<n; j++) {
                s[j] = x[j]-xprev[j];
                y[j] = g[j]-gprev[j];
            }
            invH.add_correction(s, y);

            if ((fprev-f)/std::max(std::max(std::abs(fprev), std::abs(f)), 1.)<=ftol) break;
            double gmax = 0.;
#pragma omp parallel for reduction(max:gmax)
            for (double gi : g)
                gmax = std::max(gmax, std::abs(gi));
            if (gmax <= gtol) break;
        }
    }

}

