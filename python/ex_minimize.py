#!/usr/bin/env python
"""ex_minimize.py

This examples shows how to use the C++ optimize package to find the minimum of
a python function via the pybindings interface.

"""

import time
import math
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
    avg_time = [0., 0.]
    std_time = [0., 0.]
    avg_dist = [0., 0.] # Euclidean distance between estimate and truth
    std_dist = [0., 0.]
    for iter in range(niter):
        if niter <= 10 or iter % (niter / 10) == 0:
            print 'Starting iteration',iter,'/',niter
        start = time.time()
        res1 = cppopt.minimize(myf, x0)
        end = time.time()
        delta = end-start
        dist  = math.sqrt(math.pow(res1[0] - 0.5, 2) + \
                          math.pow(res1[1] - 2., 2) + \
                          math.pow(res1[2] + 1, 2))
        avg_time[0] += delta
        avg_dist[0] += dist
        if iter:
            std_time[0] += math.pow(iter * delta - avg_time[0], 2)/ \
                           (iter * (iter  + 1))
            std_dist[0] += math.pow(iter * dist - avg_dist[0], 2)/ \
                           (iter * (iter  + 1))
        start = time.time()
        res2 = sciopt.fmin(myf, x0, disp=False)
        end = time.time()
        delta = end-start
        dist  = math.sqrt(math.pow(res2[0] - 0.5, 2) + \
                          math.pow(res2[1] - 2., 2) + \
                          math.pow(res2[2] + 1, 2))
        avg_time[1] += delta
        avg_dist[1] += dist
        if iter:
            std_time[1] += math.pow(iter * delta - avg_time[1], 2)/ \
                           (iter * (iter  + 1))
            std_dist[1] += math.pow(iter * dist - avg_dist[1], 2)/ \
                           (iter * (iter  + 1))
    print 'number of iterations =',niter
    std = math.sqrt(std_time[0])/niter
    print 'C++   avg time = %7.6f +- %7.6f' % (avg_time[0]/niter,std)
    std = math.sqrt(std_dist[0])/niter
    print 'C++   avg dist = %7.6f +- %7.6f' % (avg_dist[0]/niter,std)
    std = math.sqrt(std_time[1])/niter
    print 'scipy avg time = %7.6f +- %7.6f' % (avg_time[1]/niter,std)
    std = math.sqrt(std_dist[1])/niter
    print 'scipy avg dist = %7.6f +- %7.6f' % (avg_dist[1]/niter,std)

if __name__ == "__main__":
    main()
