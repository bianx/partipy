#!/usr/bin/env python3

import scipy.integrate
import numpy
import math

n = 100
eps = 0.5
M = 10
L = 1
t1 = 4.0

def f(z, t):
    cos = math.cos
    sin = math.sin
    cosh = math.cosh
    sinh = math.sinh
    pi = math.pi
    x = z[:n]
    y = z[n:]
    fx = numpy.zeros(n)
    fy = numpy.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            den = cosh(2*pi*dy) - cos(2*pi*dx) + eps*eps
            fx[i] -= sinh(2*pi*dy)/den
            fy[i] += sin(2*pi*dx)/den
    return numpy.hstack((fx, fy))/(2.0*n)

x0 = numpy.linspace(0, L, n)
x = x0 + 0.01 * numpy.sin(2 * math.pi * x0)
y = -0.01 * numpy.sin(2 * math.pi * x0)
z0 = numpy.hstack((x, y))
t = numpy.linspace(0, t1, M)
z, info = scipy.integrate.odeint(f, z0, t, full_output=True)
x = z[:, :n]
y = z[:, n:]

for t in range(M):
    if t > 0:
        print("")
    for i in range(n):
        print(x[t, i], y[t, i])
