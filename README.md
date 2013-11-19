optimize
========

C++ package with a list of common optimization algorithms.
This package heavily relies on the C++ Boost library (version >= 1.41.0).
To compile, you should also have installed CMake (version >= 1.26).

Install
=======

Latest sources can be checked out via:

    git clone git://github.com:bruneli/optimize.git

or (if you have write privileges):

    git clone git@github.com:bruneli/optimize.git

Then to compile

    mkdir build
    cd build
    cmake ..
    make all test

What kind of C++ function can be minimized ?
============================================

Every function which takes as an input argument a `double *` pointer and
returns a `double` can be minimized.
For example a C-like function should be something like:

    double myfunc(double * const x)
    {
      return (x[0] - 0.5)*(x[0] - 0.5) + (x[1] - 2.)*(x[1] - 2.);
    }
    double x0[2] = {0., 0.};
    minimize_nm min_nm(&myfunc, 2, (double*)x0);
    min_nm.minimize();

Using the optimize package with Python
======================================

If you have `python` and `boost-python` properly compiled, it is possible to
run the different minimization methods with a python function via a callback
function as shown in the following example:

    >>> import liboptimize as cppopt
    >>> myf = lambda x : (x[0] - 2.)*(x[0] - 2.) + (x[1] - 2.)*(x[1] - 2.) 
    >>> x0 = (0., 0.)
    >>> cppopt.minimize(myf, x0)

The file `python/ex_minimize.py` provides a comparison in term of computation 
time between the C++ `optimize` and the python based `scipy.optimize` packages
when running the Nelder-Mead (Downhill simplex) and the BFGS algorithms.

To run with the C++ optimize package, you should make sure that the PYTHONPATH
environment variable can access the `liboptimize.so` library like:

    export PYTHONPATH=/path/to/my/package/optimize/lib:${PYTHONPATH}
