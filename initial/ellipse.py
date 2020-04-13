#!/usr/bin/env python3

import math
import sys
import scipy.special
import scipy.optimize

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

def punto(x, y, lines):
    m = len(x)
    A = math.pi * a * b
    dA = A / m
    for i in range(len(x)):
        print("%.16e %.16e %.16e" % (x[i], y[i], ksi * dA))

def skel(x, y, lines):
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

ksi = 10.61
a = 1.6
b = 0.8
m0 = 4
da = a/n
db = b/n
x = []
y = []


j = 0
lines = [ ]
for i in range(n):
    a0 = da * (i + 1/2)
    b0 = db * (i + 1/2)
    m = m0 * (2 * i + 1)
    line = [ ]
    for k in range(m):
        sol = scipy.optimize.root_scalar(s, (b0, a0, k, m), x0=0, x1=2*math.pi)
        t = sol.root
        x.append(b0*math.sin(t))
        y.append(a0*math.cos(t))
        line.append(j)
        j += 1
    lines.append(line)

write(x, y, lines)
