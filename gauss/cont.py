import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy
import math

d = numpy.loadtxt("000000.dat")
x = d[:, 0]
y = d[:, 1]
o = d[:, 2]

om = numpy.zeros_like(o)

n = len(x)
d = 0.08
for i in range(n):
    for j in range(n):
        dx = x[i] - x[j]
        dy = y[i] - y[j]
        r = math.sqrt(dx*dx + dy*dy)
        d2 = d * d
        om[i] += math.exp(-r/d2) / math.pi / d2

plt.tricontour(x, y, om, 10, linewidths=0.5, colors='k')
plt.show()
