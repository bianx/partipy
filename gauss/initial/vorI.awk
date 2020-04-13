#!/bin/sh

awk '
BEGIN {
n = 100
for (i = 0; i < n; i++) {
    r = (i + 1/2) / n
    print r, vorI(r)
}
}

function f(z) {
    q = 2.56085
    return exp(-(q/z)*exp(1/(z - 1)))
}

function vorI(r) {
    Ksi = 20
    return Ksi * (1 - f(r))
}

'
