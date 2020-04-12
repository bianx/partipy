set size sq
unset key

plot "<./info -t 0.5 -p 0.1 0.5 -o twig < data/points" w l lw 3, \
     "<./info -t 0.5 -p 0.1 0.5 -o leaf < data/points" w l lw 3, \
     "<./info -t 0.5 -p 0.1 0.5 -o force < data/points" w l lw 1, \
     "data/points" w p pt 7, \
     "<echo 0.1 0.5" w p pt 7
