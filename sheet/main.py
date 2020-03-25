#!/usr/bin/env python3

import scipy.integrate
import numpy
import math
import sys

M = 10
L = 1
t1 = 4.0
me = "sheet/main.py"

eps = n = None
while True:
    sys.argv.pop(0)
    if len(sys.argv) == 0 or len(sys.argv[0]) < 2 or sys.argv[0][0] != '-':
        break
    if sys.argv[0][1] == 'n':
        sys.argv.pop(0)
        if len(sys.argv) == 0:
            sys.stderr.write("%s: -n need an argument\n" % me)
            sys.exit(2)
        n = int(sys.argv[0])
    elif sys.argv[0][1] == 'e':
        sys.argv.pop(0)
        if len(sys.argv) == 0:
            sys.stderr.write("%s: -e need an argument\n" % me)
            sys.exit(2)
        eps = float(sys.argv[0])
    else:
        sys.stderr.write("%s: unknown option '%s'\n" % (me, sys.argv[0]))
        sys.exit(2)

if n == None:
    sys.stderr.write("%s: -n is not set\n" % me)
    sys.exit(2)
if eps == None:
    sys.stderr.write("%s: -e is not set\n" % me)
    sys.exit(2)

def func(z, t):
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
z, info = scipy.integrate.odeint(func, z0, t, full_output=True)
x = z[:, :n]
y = z[:, n:]

for t in range(M):
    if t > 0:
        print("")
    for i in range(n):
        print("%.16e %.16e" % (x[t, i], y[t, i]))
