#include <iostream>
#include <deque>
#include <cassert>
#include "stlbfgs.h"

typedef std::vector<double> vector;
typedef std::deque<vector> history;

inline double dot(const vector &a, const vector &b) {
    assert(a.size()==b.size());
    double dot = 0;
#pragma omp parallel for reduction(+:dot)
    for (size_t i=0; i<a.size(); i++)
        dot += a[i]*b[i];
    return dot;
}

void two_loop(const history &S, const history &Y, const vector &g, vector &p) {
    const int m = static_cast<int>(S.size());
    const int n = static_cast<int>(g.size());
    assert(Y.size() == m);

    p = g;//vector(n);
    /*
#pragma omp parallel for
    for (int j=0; j<n; j++)
        p[j] = g[j];
    */

    for (int i=0; i<m; i++) {
        const vector &y = Y[i];
        const vector &s = S[i];
        assert(n==s.size() && n==y.size());
        double sp = dot(s, p);
        double sy = dot(s, y);
        assert(std::abs(sy)>1e-14);

#pragma omp parallel for
        for (int j=0; j<n; j++)
            p[j] -= (sp*y[j])/sy;
    }

    if (m>0) {
        double yy = dot(Y[0], Y[0]);
        assert(std::abs(yy)>1e-14);
        double c = dot(S[0], Y[0])/yy;
#pragma omp parallel for
        for (int j=0; j<n; j++)
            p[j] *= c;
    }
//    for (double v : y) std::cerr << v << std::endl;


    for (int i=m; i--; ) {
        const vector &y = Y[i];
        const vector &s = S[i];
        assert(n==s.size() && n==y.size());
        double sy = dot(s, y);
        double a = dot(s, p)/sy;
        double b = dot(y, p)/sy;
#pragma omp parallel for
        for (int j=0; j<n; j++)
            p[j] += (a-b)*s[j];
    }
}

void LBFGSopt::go(vector& x) {
    const int n = static_cast<int>(x.size());

    vector g(n), p(n), xprev(n), gprev(n);
    double f;

    history S, Y;
    for (int i=0; i<max_iter_; i++) {
        func_grad_(x, f, g);
        two_loop(S, Y, g, p);
        gprev = g;
        xprev = x;

//    for (double v : g)
//        std::cerr << v << std::endl;

    for (double v : p) std::cerr << v << std::endl;

        for (int j=0; j<n; j++) {
            x[j] += .1 * p[j];
        }

        func_grad_(x, f, g);

        S.emplace_front(n);
        Y.emplace_front(n);
        vector &s = S.front();
        vector &y = Y.front();
#pragma omp parallel for
        for (int j=0; j<n; j++) {
            s[j] = x[j]-xprev[j];
            y[j] = g[j]-gprev[j];
        }
        if (S.size()>history_depth_) S.pop_back();
        if (Y.size()>history_depth_) Y.pop_back();

        std::cerr << "f: " << f << std::endl;
        if (i==1) break;
    }
}

