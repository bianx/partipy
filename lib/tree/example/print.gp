#!/bin/sh

me=example/print.gp

for c in 1 2 4 8 16 32 64
do
    f=`printf print.%02d.png $c`
    printf %s\\n "$me: $f"
gnuplot <<!
set term pngcairo
set output "$f"
set size sq
unset key
unset border
unset xtics
unset ytics
set xrange [-1:1]
set yrange [-1:1]

plot \
     "<./print -c $c < data/points" w l lw 2, \
     "data/points" w p pt 7 ps 1.5
!
done
