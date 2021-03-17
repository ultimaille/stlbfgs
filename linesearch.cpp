#include <iostream>
#include <functional>
//#include <vector>
#include <cmath>
#include <limits>
#undef NDEBUG
#include <cassert>

typedef std::function<void(const double alpha, double &f, double &g)> func_deriv_eval;
template <typename T> auto square(const T &number) { return number * number; }

// f0 Value of function at starting point of line search.
// d0 Directional derivative at starting point of line search.
// alpha the step length.
// fa Value of function at alpha.
// mu the sufficient decrease constant. Should take a value between 0 and 1.
bool sufficient_decrease(double f0, double d0, double alpha, double fa, double mu) {
    return fa <= f0 + mu*alpha*d0;
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

double trial_value(double al, double fl, double gl, double at, double ft, double gt, double au, double fu, double gu) {
    assert(al<au); assert(al<at); assert(at<au); // TODO do we have this???
    double ac = find_cubic_minimizer(al, fl, gl, at, ft, gt);
    double aq = find_quadratic_minimizer(al, fl, gl, at, ft);
    if (ft > fl) // Case 1: a higher function value. The minimum is bracketed.
        return (std::abs(ac - al) < std::abs(aq - al)) ? ac : (aq + ac)/2.;

    double as = find_quadratic_minimizer(al, gl, at, gt);
    if (gt*gl < 0) // Case 2: A lower function value and derivatives of opposite sign. The minimum is bracketed.
        return (std::abs(ac - at) >= std::abs(as - at)) ? ac : as;

    constexpr double delta = .66; // the magic constaint is used in the Moré-Thuente paper without explanation (Section 4, Case 3).
    if (std::abs(gt) < std::abs(gl)) { // Case 3: A lower function value, derivatives of the same sign, and the magnitude of the derivative decreases.
        const double res = (std::abs(ac - at) < std::abs(as - at)) ? ac : as;
        return at > al ?
            std::min(at + delta*(au - at), res) :
            std::max(at + delta*(au - at), res);
    }

    // Case 4: A lower function value, derivatives of the same sign, and the magnitude of the derivative does not decrease.
    double ae = find_cubic_minimizer(at, ft, gt, au, fu, gu);
    return ae;
}

void line_search(const func_deriv_eval phi, const double a0, const double mu, const double eta) {
    double f0, g0;
    phi(a0, f0, g0);

    double a_l = 0.;
    double a_u = 0.;
    bool brackt = false; // set to true when a minimizer has been bracketed in an interval of uncertainty  with endpoints stx and sty
    bool stage1 = true;  // use function psi instead if phi

    // TODO alpha_min alpha_max nfev
    while (1) {
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

