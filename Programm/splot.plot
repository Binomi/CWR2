reset
set terminal epslatex color
set output 'salto1rk4.tex'
set xlabel '$l_2$ [m]'
set ylabel '$m_2$ [kg]'

set xrange [0:10]
set yrange [0:10]

set pm3d map
splot 'scan_l2m2.txt' u 1:2:5 notitle

set output
!epstopdf salto1rk4.eps

set output 'salto2rk4.tex'
splot 'scan_l2m2.txt' u 1:2:6 notitle
set output
!epstopdf salto2rk4.eps

set output 'energie_rk4.tex'
splot 'scan_l2m2.txt' u 1:2:11 notitle
set output
!epstopdf energie_rk4.eps

set output 'energie0_rk4.tex'
splot 'scan_l2m2.txt' u 1:2:9 notitle
set output
!epstopdf energie0_rk4.eps