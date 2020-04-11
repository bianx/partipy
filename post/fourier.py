#!/usr/bin/env python3

import math
import cmath
import numpy
import numpy.fft

D = numpy.loadtxt("q")
x = D[:, 0]
y = D[:, 1]
ksi = D[:, 3]
ksi0 = D[:, 4]


N = len(x)
G = numpy.arange(N) / N
p = x + 1j * y - G
f = numpy.fft.fft(p)
g = numpy.linspace(0, max(G), 2*len(G))

for n in [G[0], G[1], G[2]]:
    s = 0
    for k in range(N):
        s += 1/N * f[k] * cmath.exp(1j * 2 * math.pi * k * n)
    s += n
    print(s.real, s.imag)
