#!/bin/env python3

import numpy
import sys
import math

def triginterp(xi,x,y):
    N = len(x)
    h = 2 / N
    scale = (x[1] - x[0]) / h
    x  = x / scale
    xi = xi / scale
    P = numpy.zeros_like(xi)
    for k in range(N):
        P = P + y[k]*trigcardinal(xi-x[k], N)
    return P

def trigcardinal(x, N):
  if N % 2 == 1:
    tau = numpy.sin(N*math.pi*x/2) / (N*numpy.sin(math.pi*x/2))
  else:
    tau = numpy.sin(N*math.pi*x/2) / (N*numpy.tan(math.pi*x/2))
  tau[x==0] = 1
  return tau

D = numpy.loadtxt(sys.stdin)
x = D[:, 0]
y = D[:, 1]

N = len(x)
p = x + 1j * y

G = numpy.linspace(0, N - 1, N)
g = numpy.linspace(0, N - 1, 100 * N)
p0 = triginterp(g, G, p)

x0 = numpy.real(p0)
y0 = numpy.imag(p0)

for i in range(len(x0)):
    print(x0[i], y0[i])
