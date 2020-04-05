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

#def p(y, x):
#    r2 = x*x + y*y
#    d2 = d * d
#    return math.exp(-r2/d2) / math.pi / d2

def J2(z):
    return scipy.special.jv(2, z)
def p(y, x):
    eps = 1e-3
    d2 = d * d
    r = math.sqrt(x**2 + y**2)
    ans = (4*J2(2*r) - J2(r))/r**2 if r > eps else 15/8 - 21*r**2/32
    return ans / (3 * math.pi) / d2

def ellipse(y, x, a, b):
    return 1 if (x/b)**2 + (y/a)**2 < 1 else 0

d = 0.4
Ksi = 10.61
a = 3
b = 4
n = 7
m0 = 4
da = a/n
db = b/n
x = []
y = []
for i in range(n):
    a0 = da * (i + 1/2)
    b0 = db * (i + 1/2)
    m = m0 * (2 * i + 1)
    for k in range(m):
        sol = scipy.optimize.root_scalar(s, (a0, b0, k, m), x0=0, x1=2*math.pi)
        t = sol.root
        x.append(b0*math.sin(t))
        y.append(a0*math.cos(t))

epsabs = 0
epsrel = 1e-8
area = math.pi * a * b
sys.stderr.write("area: %g\n" % area)
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
