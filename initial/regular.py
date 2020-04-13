#!/usr/bin/env python3

import math
import sys
import scipy.special
import scipy.optimize
import scipy.integrate

def usg():
    sys.stderr.write("""\
%s [-n int] [-o skel|punto] > file
""" % me)
    sys.exit(2)

def arc(t, a, b):
    m = 1 - b**2/a**2
    return a * scipy.special.ellipeinc(t, m)

def s(t, a, b, k, n):
    c = arc(2*math.pi, a, b)
    return arc(t, a, b) - k * c / n

def punto(x, y, ksi, lines):
    m = len(x)
    for i in range(len(x)):
        print("%.16e %.16e %.16e" % (x[i], y[i], ksi[i]))

def f(z, q):
    if z == 0:
        return 0
    else:
        return math.exp(-(q/z)*math.exp(1/(z - 1)))
def vorI(r):
    if r < 1:
        return Ksi * (1 - f(r, q))
    else:
        return 0
def vorII(r):
    if r < 1:
        return Ksi * (1 - (r)**4)
    else:
        return 0
def const(r):
    if r < 1:
        return Ksi
    else:
        return 0

def skel(x, y, ksi, lines):
    nv = len(x)
    np = len(lines)
    print("SKEL")
    print("%d %d" %  (nv, np))
    for i in range(nv):
        print("%.16e %.16e %.16e" % (x[i], y[i], 0))
    for i in range(np):
        line = lines[i]
        nv = len(line) + 1
        print("%d" % nv, end = '')
        for v in line:
            print(" %d" % v, end = '')
        print(" %d" % line[0])

me = "initial/ellipse.py"
n = 8
Write = { "punto" : punto, "skel" : skel }
write = Write["punto"]
while True:
    sys.argv.pop(0)
    if len(sys.argv) == 0 or len(sys.argv[0]) < 2 or sys.argv[0][0] != '-':
        break
    if sys.argv[0][1] == 'h':
        usg()
    elif sys.argv[0][1] == 'n':
        sys.argv.pop(0)
        if len(sys.argv) == 0:
            sys.stderr.write("%s: -n needs an argument\n" % me)
            sys.exit(1)
        n = int(sys.argv[0])
    elif sys.argv[0][1] == 'o':
        sys.argv.pop(0)
        if len(sys.argv) == 0:
            sys.stderr.write("%s: -o needs an argument\n" % me)
            sys.exit(1)
        if sys.argv[0] in Write:
            write = Write[sys.argv[0]]
        else:
            sys.stderr.write("%s: unknown output type '%s'\n" % (me, sys.argv[0]))
            sys.exit(1)
    else:
        sys.stderr.write("%s: wrong option '%s'\n" % (me, sys.argv[0]))
        sys.exit(1)
Ksi = 20
q = 2.56085
a = 1.6
b = 0.8
m0 = 4
dr = 1/n
epsabs = 0
epsrel = 1e-3;

j = 0
lines = [ ]
x = []
y = []
ksi = []

fun = lambda r : r * vorI(r)
ks, err = scipy.integrate.quad(fun, 0, 1)
ks *= 2 * math.pi
sys.stderr.write("total: %g\n" % (ks * a * b));

for i in range(n):
    rlo = dr * i
    rhi = dr * (i + 1)
    a0 = dr * (i + 1/2) * a
    b0 = dr * (i + 1/2) * b
    m = m0 * (2 * i + 1)
    line = [ ]
    for k in range(m):
        tlo = scipy.optimize.root_scalar(s, (b0, a0, k, m), x0=0, x1=2*math.pi).root
        thi = scipy.optimize.root_scalar(s, (b0, a0, k + 1, m), x0=0, x1=2*math.pi).root
        
        fun = lambda r : r * vorI(r)
        fun_r = lambda r : r * r * vorI(r)
        
        r0, err = scipy.integrate.quad(fun_r, rlo, rhi)
        p0 = (thi + tlo) / 2
        ks, err = scipy.integrate.quad(fun, rlo, rhi)
        r0 /= ks
        
        x.append(b * r0 * math.sin(p0))
        y.append(a * r0 * math.cos(p0))
        ksi.append((thi - tlo) * a * b * ks)
        line.append(j)
        j += 1
    lines.append(line)

write(x, y, ksi, lines)
