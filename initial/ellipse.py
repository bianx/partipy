#!/usr/bin/env python3

import math
import sys
import scipy.special
import scipy.optimize

def arc(t, a, b):
    m = 1 - b**2/a**2
    return a * scipy.special.ellipeinc(t, m)

def s(t, a, b, k, n):
    c = arc(2*math.pi, a, b)
    return arc(t, a, b) - k * c / n

ksi = 10.61
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
        x.append(a0*math.sin(t))
        y.append(b0*math.cos(t))

m = len(x)
A = math.pi * a * b
dA = A / m
for i in range(len(x)):
    print("%.16e %.16e %.16e" % (x[i], y[i], ksi * dA))
