#!/bin/sh

me=example/box.gp

p="-0.1 -0.1"
for c in 1 2 4 8 16 32 64
do
    f=`printf box.%02d.png $c`
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
     "<echo $p" w p pt 7 ps 1.5, \
     "<./box -c $c -p $p< data/points" w l lw 4
!
done
