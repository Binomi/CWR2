reset
set terminal epslatex color
set output 'geschw1_xvx1.tex'

set xlabel '$x_1$ [m]'
set ylabel '$\dot{x_1}$ [m/s]'
set key outside bot center box

p 'rk2rk4_1.txt' u 10:14 t 'Runge-Kutta 4.Ordnung', '' u 2:6 t 'Runge-Kutta 2.Ordnung' 

set output
!epstopdf geschw1_xvx1.eps

set output 'geschw1_xvy1.tex'
set xlabel '$x_1$ [m]'
set ylabel '$\dot{y_1}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 10:15 t 'Runge-Kutta 4.Ordnung', '' u 2:7 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_xvy1.eps

set output 'geschw1_yvx1.tex'
set xlabel '$y_1$ [m]'
set ylabel '$\dot{x_1}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 11:14 t 'Runge-Kutta 4.Ordnung', '' u 3:6 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_yvx1.eps

set output 'geschw1_yvy1.tex'
set xlabel '$y_1$ [m]'
set ylabel '$\dot{y_1}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 11:15 t 'Runge-Kutta 4.Ordnung', '' u 3:7 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_yvy1.eps