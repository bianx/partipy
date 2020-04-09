import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy
import math

d = numpy.loadtxt("000020.dat")
x = d[:, 0]
y = d[:, 1]
o = d[:, 2]

plt.tricontour(x, y, o, levels = [0.5, 1, 10], linewidths=0.5, colors='k')
plt.show()

