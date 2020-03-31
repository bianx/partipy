#!/usr/bin/env python3

import scipy.integrate

def indicator(x, y):
    return 1.0 if x ** 2 + y ** 2 - 1 < 0 else 0.0

def xavg(x, y):
    return x * indicator(x, y)

def yavg(x, y):
    return y * indicator(x, y)

n = 5
h = 1.0 / n
for i in range(n):
    for j in range(n):
        xl = i * h
        xh = (i + 1) * h
        yl = j * h
        yh = (j + 1) * h
        ksi, err = scipy.integrate.dblquad(indicator, xl, xh, yl, yh)
        x, err = scipy.integrate.dblquad(xavg, xl, xh, yl, yh)
        y, err = scipy.integrate.dblquad(yavg, xl, xh, yl, yh)
        if ksi > 1e-8:
            x /= ksi
            y /= ksi
            print("%.16e %.16e %.16e" % (x, y, ksi))
            print("%.16e %.16e %.16e" % (x, -y, ksi))
            print("%.16e %.16e %.16e" % (-x, y, ksi))
            print("%.16e %.16e %.16e" % (-x, -y, ksi))            
