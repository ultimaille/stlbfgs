#include <tuple>
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
        if (D<=0) return std::numeric_limits<double>::max(); // no minumum
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

    std::tuple<double,int> trial_value(const Sample &l, const Sample &t, const Sample &u, const bool bracketed) {
        assert(l.a<t.a); assert(t.a<u.a); assert(l.d<0);
        std::cerr << l.d << " d " << t.d << std::endl;
        double ac = find_cubic_minimizer(l.a, l.f, l.d, t.a, t.f, t.d);
        if (t.f > l.f) { // Case 1: a higher function value. The minimum is bracketed.
            double aq = find_quadratic_minimizer(l.a, l.f, l.d, t.a, t.f);
//          std::cerr << "ac: " << ac << " aq: " << aq << std::endl;
//          std::cerr << "lf: " << l.f << " tf: " << t.f <<std::endl;
            double res = (std::abs(ac - l.a) < std::abs(aq - l.a)) ? ac : (aq + ac)/2.;
            assert(l.a<res); assert(res<t.a);
            return std::make_tuple(res, 1);
        }

        double as = find_quadratic_minimizer(l.a, l.d, t.a, t.d);
        if (t.d > 0) { // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
            double res = (std::abs(ac - t.a) >= std::abs(as - t.a)) ? ac : as;
            assert(l.a<res); assert(res<t.a);
            return std::make_tuple(res, 2);
        }

//      constexpr double delta = .66; // the magic constant is used in the Moré-Thuente paper (Section 4, Case 3).
        if (std::abs(t.d) < std::abs(l.d)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
//        std::cerr << "ac: " << ac << " as: " << as << std::endl;
//            bool inf_right = ((l.d + t.d)*(t.a - l.a) > 2.*(t.f - l.f)); // the cubic tends to infinity in the direction of the step
//            double res = (inf_right && (std::abs(ac - t.a) > std::abs(as - t.a))) ? ac : as; // TODO WHY > in O'Leary's code? It is < in the paper
            double res = bracketed ?
                (std::abs(ac - t.a) < std::abs(as - t.a) ? ac : as) :
                (std::abs(ac - t.a) > std::abs(as - t.a) ? ac : as);
//              std::cerr << res << std::endl;

//          std::cerr << inf_right << " " << ac << " " << as << " " << res << std::endl;
    //        std::cerr << (t.a ==res) << std::endl;
            assert(t.a<=res);
//          res = std::min(t.a + delta*(u.a - t.a), res);
//          assert(t.a<=res); assert(res<u.a);
            return std::make_tuple(res, 3);
        }

        // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
        double res = u.a;
        if (bracketed) {
            res = find_cubic_minimizer(t.a, t.f, t.d, u.a, u.f, u.d);
        }
        return std::make_tuple(res, 4);
    }

    bool line_search(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta) {
        int nfev = 0;
        auto psi = [&](const Sample &phia) {
            return Sample{ phia.a, phia.f-phi0.f-mu*phi0.d*phia.a,  phia.d-mu*phi0.d };
        };

        bool stage1 = true;  // use function psi instead if phi
        bool bracketed = false;

        Sample phil = phi0;
        Sample phit = phi(at);
        Sample phiu = phi0;
        ++nfev;

        // TODO alpha_min alpha_max nfev
        for (int i=0; i<50; i++) {
            if (!bracketed) {
              phiu.a = phit.a + 4.*(phit.a - phil.a);
//              phiu = phi(phit.a + 4.*(phit.a - phil.a));
//              ++nfev; // TODO remove eval here
            }

            std::cerr << "[" <<  phil.a << " " << phit.a << " " << phiu.a << "]" << std::endl;
            assert(phil.a<phit.a); assert(phit.a<phiu.a);


            if (sufficient_decrease(phi0, phit, mu) && curvature_condition(phi0, phit, eta)) {
                std::cerr << "Yahoo! " << nfev << std::endl << std::endl;
                return true;
            }

            // Pick next step size by interpolating either phi or psi depending on which update algorithm is currently being used.
            Sample psit = psi(phit);

            // Decide if we want to switch to using a "Modified Updating Algorithm" (shown after theorem 3.2 in the paper)
            // by switching from using function psi to using function phi.
            // The decision follows the logic in the paragraph right before theorem 3.3 in the paper.
            if (stage1 && psit.f<=0 && phit.d>0) {
                stage1 = false;
                std::cerr << "Stage 2\n";
            }

            int caseno;

            std::tie(at,caseno) = stage1 ?
                trial_value(    phil,  phit,     phiu,  bracketed) :
                trial_value(psi(phil), psit, psi(phiu), bracketed);
            bracketed = bracketed || (caseno==1 || caseno==2);
//          assert(phil.a <= at); assert(at <= phiu.a);
            std::cerr << "Case: " << caseno << std::endl;
            std::cerr << "Bracketed: " << bracketed << std::endl;

            double width = phiu.a - phil.a;
            if (1==caseno || 2==caseno) {
                phiu = phit;
            } else {
                phil = phit;
            }

//          std::cerr << phil.a << "-" << phiu.a << std::endl;

            at = std::max(phil.a, std::min(phiu.a, at));

            // Force a sufficient decrease in the size of the interval of uncertainty.
//          if (bracketed && (caseno==1 || caseno==3) && phiu.a-phil.a >= .66*width) {
//              at = (phil.a+phiu.a)/2.;
//              std::cerr << "Forcing\n";
//          }

            phit = phi(at);
            nfev++;
        }
        return false;
    }

}

