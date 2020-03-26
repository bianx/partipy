reset
set terminal postscript eps

n=400
eps="025"

t1=4
set xrange [0:t1]
#set yrange [-1.0e-3:0]

set xlabel 'time'
set ylabel 'relative Hamiltonian'

#for eps=0.5
#H0=-0.0001511141381445326

#for eps=0.25
H0=0.01329865322619408

set size 0.7, 0.7
set output "post_ham_eps".eps.".eps"

file_in="Ham_n".n."_eps".eps."_c.dat"

plot \
file_in u 1:(($2-H0)/H0) w lp t '(H(t)-H(0))/H(0)'

set terminal qt
replot
reset
