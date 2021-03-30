#include <tuple>
#include <iostream>
#include <cmath>
#include <limits>
#undef NDEBUG
#include <cassert>

#include "linesearch.h"

namespace STLBFGS {
    bool sufficient_decrease(Sample phi0, Sample phia, double mu) {
        return phia.f <= phi0.f + mu*phia.a*phi0.d;
    }

    bool curvature_condition(Sample phi0, Sample phia, double eta) {
        return std::abs(phia.d) <= eta*std::abs(phi0.d);
    }

    // minimizer of the cubic function that interpolates f(a), f'(a), f(b), f'(b) within the given interval
    double find_cubic_minimizer(double a, double fa, double ga, double b, double fb, double gb) {
//      std::cerr << a << " : " << b << std::endl;
//      std::cerr << fa << " : " << fb << std::endl;
//      std::cerr << ga << " : " << gb << std::endl;
        if (a>b) {
            std::swap( a,  b);
            std::swap(fa, fb);
            std::swap(ga, gb);
        }
        double z = 3.*(fa - fb)/(b - a) + ga + gb;
        double D = z*z - ga*gb;
        if (D<=0) return std::numeric_limits<double>::max(); // no minumum in the interval, +inf here because of the linesearch nature
        double w = std::sqrt(D); // this code assumes a<b; negate this value if b<a.
//      std::cerr << "res: " << b - ((b - a)*(gb + w - z))/(gb - ga + 2.*w);
        return b - ((b - a)*(gb + w - z))/(gb - ga + 2.*w);
    }

    // minimizer of the quadratic function that interpolates f(a), f'(a), f(b) within the given interval
    double find_quadratic_minimizer(double a, double fa, double ga, double b, double fb) {
        return a + (((b - a)*(b - a))*ga)/(2.*(fa - fb + (b - a)*ga));
    }

    // minimizer of the quadratic function that interpolates f'(a), f'(b) within the given interval
    // N.B. the function itself is undetermined since we miss information like f(a) or f(b); however the minimizer is well-defined
    double find_quadratic_minimizer(double a, double ga, double b, double gb) {
        return b + ((b - a)*gb)/(ga - gb);
    }

    std::tuple<double,int> trial_value(const Sample &l, const Sample &t, const Sample &u, const bool bracketed) {
        assert(
                (l.a<u.a && l.a<t.a && t.a<u.a && l.d<0) ||
                (l.a>u.a && l.a>t.a && t.a>u.a && l.d>0)
              );

        std::cerr << "a: " << l.a << " " << t.a << std::endl;
        std::cerr << "f: " << l.f << " " << t.f << std::endl;
        std::cerr << "d: " << l.d << " " << t.d << std::endl;
        double ac = find_cubic_minimizer(l.a, l.f, l.d, t.a, t.f, t.d);
        if (t.f > l.f) { // Case 1: a higher function value. The minimum is bracketed.
            double aq = find_quadratic_minimizer(l.a, l.f, l.d, t.a, t.f);
//          std::cerr << "ac: " << ac << " aq: " << aq << std::endl;
//          std::cerr << "lf: " << l.f << " tf: " << t.f <<std::endl;
            double res = (std::abs(ac - l.a) < std::abs(aq - l.a)) ? ac : (aq + ac)/2.;
//          assert(l.a<res); assert(res<t.a);
            return std::make_tuple(res, 1);
        }

        double as = find_quadratic_minimizer(l.a, l.d, t.a, t.d);
        if ((l.d>0 && t.d<0) || (l.d<0 && t.d>0)) { // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
            double res = (std::abs(ac - t.a) >= std::abs(as - t.a)) ? ac : as;
//          assert(l.a<res); assert(res<t.a);
            return std::make_tuple(res, 2);
        }

//      constexpr double delta = .66; // the magic constant is used in the Moré-Thuente paper (Section 4, Case 3).
        if (std::abs(t.d) <= std::abs(l.d)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
      std::cerr << "ac: " << ac << " as: " << as << std::endl;
//            bool inf_right = ((l.d + t.d)*(t.a - l.a) > 2.*(t.f - l.f)); // the cubic tends to infinity in the direction of the step
//            double res = (inf_right && (std::abs(ac - t.a) > std::abs(as - t.a))) ? ac : as; // TODO WHY > in O'Leary's code? It is < in the paper
            double res = bracketed ?
                (std::abs(ac - t.a) < std::abs(as - t.a) ? ac : as) :
                (std::abs(ac - t.a) > std::abs(as - t.a) ? ac : as);
//              std::cerr << res << std::endl;

//          std::cerr << inf_right << " " << ac << " " << as << " " << res << std::endl;
    //        std::cerr << (t.a ==res) << std::endl;
//          assert(t.a<=res);
//          res = std::min(t.a + delta*(u.a - t.a), res);
//          assert(t.a<=res); assert(res<u.a);
            return std::make_tuple(res, 3);
        }

        // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
        double res = bracketed ?
            find_cubic_minimizer(t.a, t.f, t.d, u.a, u.f, u.d) :
            u.a;
        return std::make_tuple(res, 4);
    }

    bool line_search(const linesearch_function phi, const Sample phi0, double &at, const double mu, const double eta) {
        int nfev = 0;
        bool stage1 = true;  // use function psi instead if phi
        bool bracketed = false;

        Sample phil = phi0;
        Sample phiu = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        };

        // TODO alpha_min alpha_max nfev
        for (int i=0; i<20; i++) {
            if (!bracketed)
                phiu.a = at + 4.*(at - phil.a);

            Sample phit = phi(at);
            nfev++;

//            std::cerr << "<" << phit.f << ">" << std::endl;

            std::cerr << "[" <<  phil.a << " " << phit.a << " " << phiu.a << "]\t";
//          std::cerr << "<" <<  phil.f << " " << phit.f << " " << phiu.f << ">\t";

            if (sufficient_decrease(phi0, phit, mu) && curvature_condition(phi0, phit, eta)) {
                std::cerr << "Yahoo! " << nfev << std::endl << std::endl;
                return true;
            }

            // Decide if we want to switch to using a "Modified Updating Algorithm" (shown after theorem 3.2 in the paper)
            // by switching from using function psi to using function phi.
            // The decision follows the logic in the paragraph right before theorem 3.3 in the paper.
            auto psi = [phi0, mu](const Sample &phia) {
                return Sample{ phia.a, phia.f /*- phi0.f*/ - mu*phi0.d*phia.a, phia.d - mu*phi0.d };
            };

            // Pick next step size by interpolating either phi or psi depending on which update algorithm is currently being used.
            stage1 = stage1 && (psi(phit).f>0 || phit.d<=0);

            if (stage1)
                std::cerr << "stage 1 ";
            else
                std::cerr << "stage 2 ";

            int caseno;
            std::tie(at, caseno) = stage1 ?
                trial_value(psi(phil), psi(phit), psi(phiu), bracketed) : // TODO NaN!
                trial_value(    phil,      phit,      phiu,  bracketed);

            std::cerr << "<" << at << ">" << std::endl;
            bracketed = bracketed || (caseno==1 || caseno==2);
//          assert(phil.a <= at); assert(at <= phiu.a);
            std::cerr << " case " << caseno << " bracketed " << bracketed << std::endl;

//          double width = std::abs(phiu.a - phil.a);
            if (1==caseno) {
                phiu = phit;
            } else if (2==caseno) {
                phiu = phil;
                phil = phit;
            } else {
                phil = phit;
            }


            if (bracketed && (caseno==1 || caseno==3)) {
                at = phil.a<phiu.a ?
                    std::min(phil.a + .66*(phiu.a-phil.a), at) :
                    std::max(phil.a + .66*(phiu.a-phil.a), at);
            }

            // Force a sufficient decrease in the size of the interval of uncertainty.
//          if (bracketed && (caseno==1 || caseno==3) && std::abs(phiu.a-phil.a) >= .66*width) {
//              at = (phil.a+phiu.a)/2.;
//              std::cerr << "Forcing\n";
//          }

            at = phil.a<phiu.a ?
                std::max(phil.a, std::min(phiu.a, at)) :
                std::min(phil.a, std::max(phiu.a, at));
          std::cerr << phil.a << " - " << phiu.a << ":" << at << std::endl;
        }
        return false;
    }

}

