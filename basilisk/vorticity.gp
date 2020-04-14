set term pngcairo enhanced size 640,640
set size ratio -1
unset key
unset xtics
unset ytics
unset border
unset colorbox
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

set output 'vorticity0.png'
set multiplot layout 2,2 scale 1.5
splot 'omega-000'
splot 'omega-002'
splot 'omega-003'
splot 'omega-008'
unset multiplot

set output 'vorticity1.png'
set multiplot layout 2,2 scale 1.5
splot 'omega-016'
splot 'omega-024'
splot 'omega-032'
splot 'omega-048'
