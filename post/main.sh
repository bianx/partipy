#!/bin/sh

for i in ../gauss/*.dat
do
    o=`basename $i .dat`.png
    ./skel ../gauss/initial/vorI.skel "$i" > /tmp/t
    gnuplot <<!
    set output "$o"
    set term pngcairo mono
    set size sq
    set key off
    unset border
    unset xtics
    unset ytics
    plot [-1.9:1.9][-1.9:1.9] "/tmp/t" w l lw 2, "$i" w p ps 0.5 pt 7
!
    if test $? -ne 0
    then exit 1
    fi
done
