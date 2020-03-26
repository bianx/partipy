# Line

## geomview

     make
     ./main -m 40 -n 400 -e 0.5 -o skel
     co.geomview -t -0.5 0 1.75 -n none -a appearence/a  -O -p cat *.skel
     convert *.ppm o.gif

## gnuplot

     make
     ./main -m 40 -n 400 -e 0.5 -o gnuplot
     ./gnuplot *.dat
     convert *.png o.gif
