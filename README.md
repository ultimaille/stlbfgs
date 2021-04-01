# Basic C++ L-BFGS implementation using plain STL

L-BFGS is a quasi-Newton optimization algorithm for solving large nonlinear optimization problems [1]. It employs function value and gradient information to search for the local optimum. It uses (as the name suggests) the BGFS (Broyden-Goldfarb-Fletcher-Shanno) algorithm to approximate the Hessian matrix. The size of the memory available to store the approximation of the Hessian is limited (hence L- in the name): approximated Hessian matrix is not stored directly, but rather built from the history of past updates. This implementation uses Moré-Thuente line search algorithm [2].

The project has zero external dependencies, no Eigen, nothing, plain standard library. It is a from-scratch implementation and not a wrapper/translation of the original Fortran/Matlab codes by Jorge Nocedal / Dianne O'Leary.

Line search routine is tested against data [2] (Tables 1-6), and L-BFGS routine is tested on problems taken from [3].

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

