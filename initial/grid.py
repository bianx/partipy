#!/usr/bin/env python3

import sys
import scipy.integrate
import math

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

def f(z, q):
    return math.exp(-(q/z)*math.exp(1/(z - 1)))

def vorI(y, x, a, b):
    r = math.sqrt((x/a)**2 + (y/b)**2)
    if r < 1:
        return Ksi * (1 - f(r, q))
    else:
        return 0

def vorII(y, x, a, b):
    r = math.sqrt((x/a)**2 + (y/b)**2)
    if r < 1:
        return Ksi * (1 - (r)**4)
    else:
        return 0

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
    a = 1
    b = 1
    nx = 5
    ny = 5
    vorticity = circle
elif shape == "vorI":
    a = 0.8
    b = 2 * 0.8
    nx = 20
    ny = 20
    vorticity = lambda x, y : vorI(x, y, a, b)
    Ksi = 20
    q = 2.56085
elif shape == "vorII":
    a = 0.8
    b = 2 * 0.8
    nx = 20
    ny = 20
    vorticity = lambda x, y : vorII(x, y, a, b)
    Ksi = 20
    q = 2.56085    
else:
    sys.stderr.write("%s: unknown shape '%s'\n" % (me, shape))
    sys.exit(2)

# h = max(a, b) / nx
dx = a / nx
dy = b / ny
for i in range(nx):
    for j in range(ny):
        xl = i * dx
        xh = (i + 1) * dx
        yl = j * dy
        yh = (j + 1) * dy
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
