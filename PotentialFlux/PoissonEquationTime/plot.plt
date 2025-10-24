#!/usr/bin/gnuplot
set terminal pngcairo size 1000, 900 enhanced dashed font 'Verdana,26'

set mxtics 2
set mytics 2
set border lw 3

set style function linespoints
set style line 1 lw 4 lc rgb 'black' ps 3 pt 6  pi 50
set style line 2 lw 4 lc rgb 'red' ps 3 pt 9  pi 50
set style line 3 lw 4 lc rgb 'blue' ps 3 pt 12 pi 50
set style line 4 lw 4 lc rgb 'green' ps 3 pt 7  pi 50


set output 'Density.png'
set xlabel 'x'
set ylabel 'Charge Density'
set key top

set xrange [0:0.05]

plot './Density_D_1e-6.txt' using 1:2 with lines ls 1 title 'Slow Diffusion',\
     './Density_D_5e-6.txt' using 1:2 with lines ls 2 title 'Fast Diffusion',\
     './Initial_rho.txt'    using 1:2 with lines ls 3 title 'Initial Rho',\
     './Density_M_High.txt'    using 1:2 with lines ls 4 title 'High Mobility',\

