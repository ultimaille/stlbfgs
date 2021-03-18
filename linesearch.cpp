#include <iostream>
#include <functional>
//#include <vector>
#include <cmath>
#include <limits>
#undef NDEBUG
#include <cassert>

typedef std::function<void(const double alpha, double &f, double &g)> func_deriv_eval;
template <typename T> auto square(const T &number) { return number * number; }

struct Sample {
    double a, f, d;
};


// f0 Value of function at starting point of line search.
// d0 Directional derivative at starting point of line search.
// alpha the step length.
// fa Value of function at alpha.
// mu the sufficient decrease constant. Should take a value between 0 and 1.
bool sufficient_decrease(Sample phi0, Sample phia, double mu) {
    return phia.f <= phi0.f + mu*phia.a*phi0.d;
}

bool curvature_condition(Sample phi0, Sample phia, double eta) {
    return std::abs(phia.d) <= eta*std::abs(phi0.d);
}

// minimizer of the cubic function that interpolates f(a), f'(a), f(b), f'(b) within the given interval ]a, b[
double find_cubic_minimizer(double a, double fa, double ga, double b, double fb, double gb) {
//  if (std::abs(a - b) < precision) return a;
    const double z = 3.*(fa - fb)/(b - a) + ga + gb;
    const double w = std::sqrt(z*z - ga*gb);
    return b - ((b - a)*(gb + w - z))/(gb - ga + 2.*w);
}

// minimizer of the quadratic function that interpolates f(a), f'(a), f(b) within the given interval ]a, b[
double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb) {
//  if (std::abs(a - b) < precision) return a;
    return a + (square(b - a)*ga)/(2.*(fa - fb + (b - a)*ga));
}

// minimizer of the quadratic function that interpolates f'(a), f'(b) within the given interval ]a, b[
// N.B. the function itself is undetermined since we miss information like f(a) or f(b); however the minimizer is well-defined
double find_quadratic_minimizer(double a, double ga, double b, double gb) {
//  if (std::abs(a - b) < precision) return a;
    return b + ((b - a)*gb)/(ga - gb);
}

double trial_value(const Sample &l, const Sample &t, const Sample &u) {
    assert(l.a<u.a); assert(l.a<t.a);// assert(t.a<u.a); // TODO do we have this???
    double ac = find_cubic_minimizer(l.a, l.f, l.d, t.a, t.f, t.d);
    double aq = find_quadratic_minimizer(l.a, l.a, l.d, t.a, t.f);
    if (t.f > l.f) // Case 1: a higher function value. The minimum is bracketed.
        return (std::abs(ac - l.a) < std::abs(aq - l.a)) ? ac : (aq + ac)/2.;

    double as = find_quadratic_minimizer(l.a, l.d, t.a, t.d);
    if (t.d*l.d < 0) // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
        return (std::abs(ac - t.a) >= std::abs(as - t.a)) ? ac : as;

    constexpr double delta = .66; // the magic constaint is used in the Moré-Thuente paper without explanation (Section 4, Case 3).
    if (std::abs(t.d) < std::abs(l.d)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
        const double res = (std::abs(ac - t.a) < std::abs(as - t.a)) ? ac : as;
        return t.a > l.a ?
            std::min(t.a + delta*(u.a - t.a), res) :
            std::max(t.a + delta*(u.a - t.a), res);
    }

    // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
    double ae = find_cubic_minimizer(t.a, t.f, t.d, u.a, u.f, u.d);
    return ae;
}

void line_search(const func_deriv_eval func, double at, double mu, double eta) {
    auto phi = [&](const double a) {
        double f, d;
        func(a, f, d);
        return Sample{a, f, d};

    };
    Sample phi0 = phi(0);

    auto psi = [&](const Sample &phia) {
        return Sample{ phia.a, phia.f-phi0.f-mu*phi0.d*phia.a,  phia.d-mu*phi0.d };
    };

    bool stage1 = true;  // use function psi instead if phi

    Sample phil = phi0;
    Sample phit = phi(at);
    Sample phiu = phi(phit.a + 4.*(phit.a - phil.a));

    // TODO alpha_min alpha_max nfev
for (int i=0; i<10; i++) {
        if (phil.a > phiu.a)
            std::swap(phil, phiu);

        std::cerr << phil.a << " " << phit.a << " " << phiu.a << std::endl;

//      assert(phil.a<phiu.a); assert(phil.a<phit.a); assert(phit.a<phiu.a); // TODO do we have this???

        if (sufficient_decrease(phi0, phit, mu) && curvature_condition(phi0, phit, eta)) {
            std::cerr << "Yahoo!" << std::endl;
            return;
        }

        // Pick next step size by interpolating either phi or psi depending on which update algorithm is currently being used.
        Sample psit = psi(phit);
        double at = stage1 ? trial_value(phil, phit, phiu) : trial_value(psi(phil), psit, psi(phiu));
        phit = phi(at);
        psit = psi(phit);

        // Decide if we want to switch to using a "Modified Updating Algorithm" (shown after theorem 3.2 in the paper)
        // by switching from using function psi to using function phi.
        // The decision follows the logic in the paragraph right before theorem 3.3 in the paper.
        if (stage1 && psit.f<=0 && phit.d>0) {
            stage1 = false;
        }

        if (stage1) {
            Sample psil = psi(phil), psiu = psi(phiu);
            if (psit.f > psil.f) {
                phiu = phit;
            } else if (psit.d*(psit.a - psil.a) > 0) {
                phil = phit;
            } else {
                phiu = phil;
                phil = phit;
            }
        } else {
            // Update the interval that will be used to generate the next step using the "Modified Updating Algorithm" (right after theorem 3.2 in the paper).
            if (phit.f > phil.f) {
                phiu = phit;
            } else if (phit.d*(phit.a - phil.a) > 0) {
                phil = phit;
            } else {
                phiu = phil;
                phil = phit;
            }
        }
    }

}




/*
typedef std::vector<double> vector;
typedef std::function<void(std::vector<double>& x, double &f, std::vector<double>& g)> func_grad_eval;


void line_search(func_grad_eval phi,
        double step0,
        double alpha = 1.,
        double c1 = 1e-4,
        double c2 = 1e-1,
        double xtol = 1e-14,
        double alpha_min = 0.,
        double alpha_max = std::numeric_limits<double>::max(),
        int maxfev = std::numeric_limits<int>::max(),
        double delta = 0.66) {

        constexpr double xtrapf = 4.;
        if (alpha <= 0) {
            std::cerr << "alpha < 0.0" << std::endl;
            return;
        }
        if (c1 < 0) {
            std::cerr << "c1 < 0.0" << std::endl;
            return;
        }
        if (c2 < 0) {
            std::cerr << "c2 < 0.0" << std::endl;
            return;
        }
        if (xtol < 0) {
            std::cerr << "xtol < 0.0" << std::endl;
            return;
        }
        if (alpha_min < 0) {
            std::cerr << "alpha_min < 0.0" << std::endl;
            return;
        }
        if (alpha_max < 0) {
            std::cerr << "alpha_max < 0.0" << std::endl;
            return;
        }
        if (maxfev <= 0) {
            std::cerr << "maxfev < 0" << std::endl;
            return;
        }
}
*/

