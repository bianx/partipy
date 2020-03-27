reset
set terminal postscript eps

n=400
eps="05"

t1=4
set xrange [0:t1]
#set yrange [-1.0e-3:0]

set xlabel 'time'
set ylabel 'total linear momenta'

#for eps=0.5

#for eps=0.25

set size 0.7, 0.7
set output "post_mom_eps".eps.".eps"

file_in="post_n".n."_eps".eps."_c.dat"

plot \
file_in u 1:3 w lp t 'Mom_x',\
file_in u 1:4 w lp t 'Mom_y'

set terminal qt
replot
reset
