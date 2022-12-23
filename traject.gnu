set term gif size 512,512 animate delay 3
set output "trajectory.gif"

set key off
set size square
unset border
unset xtics
unset ytics
unset ztics

do for [j=0:99] { 

set view equal xyz
set xrange [-L/2:L/2]
set yrange [-L/2:L/2]
set zrange [-L/2:L/2]

splot "trajectory.xyz" i j u 2:3:4:(1) w p pt 7

}