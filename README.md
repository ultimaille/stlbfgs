# Basic C++ L-BFGS implementation using plain STL

[![Build](https://github.com/ssloy/stlbfgs/actions/workflows/continuous.yml/badge.svg)](https://github.com/ssloy/stlbfgs/actions/workflows/continuous.yml) [![Nightly](https://github.com/ssloy/stlbfgs/actions/workflows/nightly.yml/badge.svg)](https://github.com/ssloy/stlbfgs/actions/workflows/nightly.yml)

L-BFGS is a quasi-Newton optimization algorithm for solving large nonlinear optimization problems [1]. It employs function value and gradient information to search for the local optimum. It uses (as the name suggests) the BGFS (Broyden-Goldfarb-Fletcher-Shanno) algorithm to approximate the Hessian matrix. The size of the memory available to store the approximation of the Hessian is limited (hence L- in the name): approximated Hessian matrix is not stored directly, but rather built from the history of past updates. This implementation uses Moré-Thuente line search algorithm [2].

The project has zero external dependencies, no Eigen, nothing, plain standard library. It is a from-scratch implementation and not a wrapper/translation of the original Fortran/Matlab codes by Jorge Nocedal / Dianne O'Leary.

Line search routine is tested against data [2] (Tables 1-6), and L-BFGS routine is tested on problems taken from [3].

# Hello world
Here is a minimal usage example: it shows a minimization of a 2D quadratic function f(x, y) = (x-7)^2 + (y-1)^2, whose minimizer is obvious (x, y) = (7, 1):

```cpp
#include "stlbfgs.h"

int main() {
    std::vector<double> x = {10, 10};

    const STLBFGS::Optimizer::func_grad_eval func = [](const std::vector<double>& x, double& f, std::vector<double>& g) {
        f = (x[0] - 7)*(x[0] - 7) +
            (x[1] - 1)*(x[1] - 1);
        g[0] = 2*(x[0] - 7);
        g[1] = 2*(x[1] - 1);
    };

    STLBFGS::Optimizer opt(2, func);
    opt.run(x);

    return std::abs(x[0]-7)>1e-3 || std::abs(x[1]-1)>1e-3x;
}

```

# Compile and run unit tests:
```sh
git clone https://github.com/ssloy/stlbfgs.git &&
cd stlbfgs &&
mkdir build &&
cd build &&
cmake -DUNIT_TESTS:BOOL=ON .. &&
make -j &&
cd tests &&
ctest
```

# References
**[1]** Dong C. Liu, Jorge Nocedal. On the limited memory BFGS method for large scale optimization. Mathematical Programming (1989).

**[2]** Jorge J. Moré, David J. Thuente. Line search algorithms with guaranteed sufficient decrease. ACM Transactions on Mathematical Software  (1994)

**[3]** Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom, "Testing Unconstrained Optimization Software", ACM Transactions on Mathematical Software  (1981)

