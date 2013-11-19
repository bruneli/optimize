#!/usr/bin/env python
"""ex_minimize.py

This examples shows how to use the C++ optimize package to find the minimum of
a python function via the pybindings interface.

"""

import math
import random
import time
import scipy.optimize as sciopt
try:
    import liboptimize as cppopt
except ImportError:
    print r"""liboptimize is not found:
- did you include optimize/lib path to your PYTHONPATH environment variable ?
- have you checked if python and boost-python are properly compiled ?
"""

def myf(x):
    val  = 1.  * math.pow(x[0] - 0.5, 2)
    val += 10. * math.pow(x[1] - 2., 2)
    val += 4.  * math.pow(x[2] + 1., 2)
    return val

def main():
    x0 = (2., 2., 2.)
    niter = 2000
    nalgos = 4
    avg_time = [0.]*nalgos
    std_time = [0.]*nalgos
    avg_dist = [0.]*nalgos # Euclidean distance between estimate and truth
    std_dist = [0.]*nalgos
    for iter in range(niter):
        if niter <= 10 or iter % (niter / 10) == 0:
            print 'Starting iteration',iter,'/',niter
        seq = range(nalgos)
        random.shuffle(seq)
        for i in seq:
            start = time.time()
            if i == 0:
                res = cppopt.minimize(myf, x0)
            elif i == 1:
                res = sciopt.fmin(myf, x0, disp=False)
            elif i == 2:
                res = cppopt.minimize(myf, x0, method="BFGS")
            else:
                res = sciopt.fmin_bfgs(myf, x0, disp=False)
            end = time.time()
            delta = end-start
            dist  = math.sqrt(math.pow(res[0] - 0.5, 2) + \
                              math.pow(res[1] - 2., 2) + \
                              math.pow(res[2] + 1, 2))
            avg_time[i] += delta
            avg_dist[i] += dist
            if iter:
                std_time[i] += math.pow(iter * delta - avg_time[i], 2)/ \
                                        (iter * (iter  + 1))
                std_dist[i] += math.pow(iter * dist - avg_dist[i], 2)/ \
                                        (iter * (iter  + 1))
    print 'number of iterations =',niter
    for i,name in enumerate(['Nelder-Mead (C++)','Nelder-Mead (scipy)',
                             'BFGS (C++)','BFGS (scipy)']):
        std = math.sqrt(std_time[i])/niter
        print '%20s avg time = %7.6f +- %7.6f' % (name,avg_time[i]/niter,std)
        std = math.sqrt(std_dist[i])/niter
        print '%20s avg dist = %7.6f +- %7.6f' % (name,avg_dist[i]/niter,std)

if __name__ == "__main__":
    main()
