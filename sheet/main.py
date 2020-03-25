#!/usr/bin/env python3

import scipy.integrate
import numpy
import math

N = 400
M = 100
L = 0.5
t1 = 4.0
eps = 0.1

def f(z, t):
    cos = math.cos
    sin = math.sin
    cosh = math.cosh
    sinh = math.sinh
    pi = math.pi
    x = z[:N]
    y = z[N:]
    fx = numpy.zeros(N)
    fy = numpy.zeros(N)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            den = cosh(2*pi*dy) - cos(2*pi*dx) + eps*eps
            fx[i] -= sinh(2*pi*dy)/den
            fy[i] += sin(2*pi*dx)/den
    return numpy.hstack((fx, fy))/(2.0*N)


x0 = numpy.linspace(0, L, N)
y0 = numpy.zeros(N)
z0 = numpy.hstack((x0, y0))
t0 = 0
t = numpy.linspace(t0, t1, M)
z, info = scipy.integrate.odeint(f, z0, t, full_output=True)
x = z[:, :N]
y = z[:, N:]

for t in range(M):
    if t > 0:
        print("")
    for i in range(N):
        print(x[t, i], y[t, i])
