set term pngcairo
set output 'level.png'
set cbrange [3:8]
set multiplot layout 2,3 scale 1.6,1.6
splot 'level-000'
splot 'level-001'
splot 'level-002'
splot 'level-003'
splot 'level-004'
splot 'level-005'
unset multiplot
