reset
set terminal epslatex color
set output 'energie4.tex'

set xlabel '$t$ [s]'
set ylabel '$E$ [J]'
set key top left

p 'rk2rk4_4.txt' u 1:19 t 'Runge-Kutta 4.Ordnung', '' u 1:18 t 'Runge-Kutta 2.Ordnung' 

set output
!epstopdf energie4.eps

reset
set terminal epslatex color
set output 'dsum4.tex'

set xlabel '$t$ [s]'
set ylabel 'dsum [m]'
set key top left

p 'rk2rk4_4.txt' u 1:20 t 'Masse 1', '' u 1:21 t 'Masse 2' 

set output
!epstopdf dsum4.eps