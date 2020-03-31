#!/usr/bin/env python3

import sys
import scipy.integrate

def circle(y, x):
    return 1 if x ** 2 + y ** 2 < 1 else 0

def cubic(y, x):
    inside = x ** 2 + y ** 2 < 1
    return (1 - (x**2 + y**2))**3 if inside else 0

def ellipse(y, x, a, b):
    return 1 if (x/a)**2 + (y/b)**2 < 1 else 0

def xavg(y, x, vorticity):
    return x * vorticity(y, x)

def yavg(y, x, vorticity):
    return y * vorticity(y, x)

me = "initial/grid.py"
shape = "ellipse"
epsabs = 0
epsrel = 1e-3
small = 1e-12

sys.argv.pop(0)
if len(sys.argv) == 0:
    sys.stderr.write("%s: shape argumen is not set\n" % me)
    sys.exit(2)
shape = sys.argv.pop(0)

if shape == "ellipse":
    a = 1
    b = 2
    nx = 5
    ny = 10
    vorticity = lambda x, y : ellipse(x, y, a, b)
elif shape == "cubic":
    nx = 5
    ny = 5
    vorticity = cubic
elif shape == "circle":
    nx = 5
    ny = 5
    vorticity = circle
else:
    sys.stderr.write("%s: unknown shape '%s'\n" % (me, shape))
    sys.exit(2)

h = 1 / nx
for i in range(nx):
    for j in range(ny):
        xl = i * h
        xh = (i + 1) * h
        yl = j * h
        yh = (j + 1) * h
        ksi, err = scipy.integrate.dblquad(vorticity, xl, xh, yl, yh, (), epsabs, epsrel)
        x, err = scipy.integrate.dblquad(xavg, xl, xh, yl, yh, (vorticity,), epsabs, epsrel)
        y, err = scipy.integrate.dblquad(yavg, xl, xh, yl, yh, (vorticity,), epsabs, epsrel)
        if ksi > small:
            x /= ksi
            y /= ksi
            print("%.16e %.16e %.16e" % (x, y, ksi))
            print("%.16e %.16e %.16e" % (x, -y, ksi))
            print("%.16e %.16e %.16e" % (-x, y, ksi))
            print("%.16e %.16e %.16e" % (-x, -y, ksi))            
