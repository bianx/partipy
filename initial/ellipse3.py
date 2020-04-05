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

def f(z, q):
    return math.exp(-(q/z)*math.exp(1/(z - 1)))

def vorI(x, y, a, b):
    r = math.sqrt((x/a)**2 + (y/b)**2)
    if r < r0:
        return Ksi * (1 - f(r/r0, q))
    else:
        return 0

def vorII(x, y, a, b):
    r = math.sqrt((x/a)**2 + (y/b)**2)
    if r < r0:
        return Ksi * (1 - (r/r0)**4)
    else:
        return 0

vor = vorI
Ksi = 20
r0 = 0.8
q = 2.56085
a = 1
b = 2
n = 8
m0 = 4
da = a/n
db = b/n
x = []
y = []
for i in range(n):
    a0 = da * (i + 1/2) * r0
    b0 = db * (i + 1/2) * r0
    m = m0 * (2 * i + 1)
    for k in range(m):
        sol = scipy.optimize.root_scalar(s, (a0, b0, k, m), x0=0, x1=2*math.pi)
        t = sol.root
        x.append(a0*math.sin(t))
        y.append(b0*math.cos(t))

m = len(x)
A = math.pi * a * b
dA = A / m
for u, v in zip(x, y):
    print(u, v, vor(u, v, a, b) * dA)
