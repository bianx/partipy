reset
set terminal postscript eps
i = 20
j = 0

t1=4
M=20
dt=4.0/20

t=4

n=400
eps="05"

round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)

i=round(t/dt)
print dt, i, t
set size 1, 0.5
set output "eps".eps."_time_".t.".eps"

L=1
L2=2
set xrange [0:L2]
set yrange [-0.275:0.275]

set xlabel 'x'
set ylabel 'y'
set xtics 1
set ytics 0.275
wi=5
si=0.5

file_in="output_n".n."_eps".eps."_c.dat"

print file_in
plot \
file_in every :::i::i u 1:2 w l lc 1 lw wi t 't='.t,\
file_in every :::i::i u ($1+L):2 w p pt 7 ps si lc 1 lw wi t ''


set terminal qt
replot
reset

#, \"output_p.dat" every :::i::i w lp t 'Python'
