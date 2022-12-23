set view equal xyz

set xrange [0:L]
set yrange [0:L]
set zrange [0:L]

splot "trajectory.xyz" i 99 u 2:3:4:(1) w p pt 7
pause -1