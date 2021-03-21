#include <iostream>
#include <cmath>
#include <limits>
#undef NDEBUG
#include <cassert>

#include "linesearch.h"

namespace STLBFGS {

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
        assert(a<b); assert(ga<0);
        double z = 3.*(fa - fb)/(b - a) + ga + gb;
        double D = z*z - ga*gb;
        if (D<0) return std::numeric_limits<double>::max(); // no minumum
        double w = std::sqrt(D);
        return b - ((b - a)*(gb + w - z))/(gb - ga + 2.*w);
    }

    // minimizer of the quadratic function that interpolates f(a), f'(a), f(b) within the given interval ]a, b[
    double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb) {
        assert(a<b); assert(ga<0);
        return a + (((b - a)*(b - a))*ga)/(2.*(fa - fb + (b - a)*ga));
    }

    // minimizer of the quadratic function that interpolates f'(a), f'(b) within the given interval ]a, b[
    // N.B. the function itself is undetermined since we miss information like f(a) or f(b); however the minimizer is well-defined
    double find_quadratic_minimizer(double a, double ga, double b, double gb) {
        assert(a<b); assert(ga<0);
        return b + ((b - a)*gb)/(ga - gb);
    }

    double trial_value(const Sample &l, const Sample &t, const Sample &u) {
        assert(l.a<u.a); assert(l.a<t.a); assert(t.a<u.a); assert(l.d<0); // TODO do we have this???
        double ac = find_cubic_minimizer(l.a, l.f, l.d, t.a, t.f, t.d);
        double aq = find_quadratic_minimizer(l.a, l.a, l.d, t.a, t.f);
        if (t.f > l.f) { // Case 1: a higher function value. The minimum is bracketed.
            double res = (std::abs(ac - l.a) < std::abs(aq - l.a)) ? ac : (aq + ac)/2.;
            assert(l.a<res); assert(res<t.a);
            return res;
        }

        double as = find_quadratic_minimizer(l.a, l.d, t.a, t.d);
        if (t.d > 0) { // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
            double res = (std::abs(ac - t.a) >= std::abs(as - t.a)) ? ac : as;
            assert(l.a<res); assert(res<t.a);
            return res;
        }

        constexpr double delta = .66; // the magic constaint is used in the Moré-Thuente paper without explanation (Section 4, Case 3).
        if (std::abs(t.d) < std::abs(l.d)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
            const double res = (std::abs(ac - t.a) < std::abs(as - t.a)) ? ac : as;
            return std::min(t.a + delta*(u.a - t.a), res);
        }

        // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
        double ae = find_cubic_minimizer(t.a, t.f, t.d, u.a, u.f, u.d);
        return ae;
    }

    void line_search(const linesearch_function phi, double at, double mu, double eta) {
        /*
           auto phi = [&](const double a) {
           double f, d;
           func(a, f, d);
           return Sample{a, f, d};
           };
         */
        Sample phi0 = phi(0);

        auto psi = [&](const Sample &phia) {
            return Sample{ phia.a, phia.f-phi0.f-mu*phi0.d*phia.a,  phia.d-mu*phi0.d };
        };

        bool stage1 = true;  // use function psi instead if phi

        Sample phil = phi0;
        Sample phit = phi(at);
        Sample phiu = phi(phit.a + 4.*(phit.a - phil.a));

        // TODO alpha_min alpha_max nfev
        for (int i=0; i<30; i++) {
            if (true) // TODO
                phiu = phi(phit.a + 4.*(phit.a - phil.a));

            std::cerr << phil.a << " " << phit.a << " " << phiu.a << std::endl;

            //      assert(phil.a<phiu.a); assert(phil.a<phit.a); assert(phit.a<phiu.a); // TODO do we have this???

            if (sufficient_decrease(phi0, phit, mu) && curvature_condition(phi0, phit, eta)) {
                std::cerr << "Yahoo!" << std::endl;
                //          return;
            }

            // Pick next step size by interpolating either phi or psi depending on which update algorithm is currently being used.
            Sample psit = psi(phit);
            double at = stage1 ?
                trial_value(    phil,  phit,     phiu) :
                trial_value(psi(phil), psit, psi(phiu));
            phit = phi(at);
            psit = psi(phit);

            // Decide if we want to switch to using a "Modified Updating Algorithm" (shown after theorem 3.2 in the paper)
            // by switching from using function psi to using function phi.
            // The decision follows the logic in the paragraph right before theorem 3.3 in the paper.
            if (stage1 && psit.f<=0 && phit.d>0) {
                stage1 = false;
            }

            Sample psil = psi(phil);
            if (
                    ( stage1 && psit.f > psil.f) ||
                    (!stage1 && phit.f > phil.f)
               ) {
                phiu = phit;
                std::cerr << 1 << std::endl;
            } else if (
                    ( stage1 && psit.d < 0) ||
                    (!stage1 && phit.d < 0)
                    ) {
                phil = phit;
                phit = phiu;
                std::cerr << 2 << std::endl;
            } else {
                phiu = phit;
                std::cerr << 3 << std::endl;
            }
        }
    }

}

