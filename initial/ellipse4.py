#!/usr/bin/env python3

import math
import sys
import scipy.special
import scipy.optimize
import scipy.integrate
import numpy

def arc(t, a, b):
    m = 1 - b**2/a**2
    return a * scipy.special.ellipeinc(t, m)

def s(t, a, b, k, n):
    c = arc(2*math.pi, a, b)
    return arc(t, a, b) - k * c / n

def p(y, x):
    r2 = x*x + y*y
    d2 = d * d
    return math.exp(-r2/d2) / math.pi / d2

def J2(z):
    return scipy.special.jv(2, z)
#def p(y, x):
#    eps = 1e-3
#    d2 = d * d
#    r = math.sqrt(x**2 + y**2)
#    ans = (4*J2(2*r) - J2(r))/r**2 if r > eps else 15/8 - 21*r**2/32
#    return ans / (3 * math.pi) / d2

def ellipse(y, x, a, b):
    return 1 if (x/b)**2 + (y/a)**2 < 1 else 0

D = numpy.loadtxt(sys.argv[1])
x = D[:, 0].tolist()
y = D[:, 1].tolist()

a = 0.8
b = 1.6
d = 0.4
epsabs = 0
epsrel = 1e-6
area = math.pi * a * b
Ksi = 20
n = len(x)
A = numpy.empty((n, n))
for i in range(n):
    sys.stderr.write("%02d of %02d\n" % (i, n))
    for j in range(i, n):
        fun = lambda v, u : p(v - y[i], u - x[i]) * p(v - y[j], u - x[j]) * ellipse(v, u, a, b)
        ans, err = scipy.integrate.dblquad(fun, -b, b, -a, a, (), epsabs, epsrel)
        A[i, j] = A[j, i] = ans

B = numpy.empty(n)
#B[n] = area * Ksi
for i in range(n):
    fun = lambda v, u : p(v - y[i], u - x[i]) * ellipse(v, u, a, b)
    ans, err = scipy.integrate.dblquad(fun, -b, b, -a, a, (), epsabs, epsrel)
    B[i] = ans * Ksi

ksi = numpy.linalg.solve(A, B)
sys.stderr.write("min, avg, max: %g %g %g\n" % (numpy.min(ksi), numpy.mean(ksi), numpy.max(ksi)))
for i in range(n):
    print("%.16e %.16e %.16e" % (x[i], y[i], ksi[i]))
