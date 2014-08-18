reset
set terminal epslatex color
set output 'geschw1_xvx2.tex'

set xlabel '$x_2$ [m]'
set ylabel '$\dot{x_2}$ [m/s]'
set key outside bot center box

p 'rk2rk4_1.txt' u 12:16 t 'Runge-Kutta 4.Ordnung', '' u 4:8 t 'Runge-Kutta 2.Ordnung' 

set output
!epstopdf geschw1_xvx2.eps

set output 'geschw1_xvy2.tex'
set xlabel '$x_2$ [m]'
set ylabel '$\dot{y_2}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 12:17 t 'Runge-Kutta 4.Ordnung', '' u 4:9 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_xvy2.eps

set output 'geschw1_yvx2.tex'
set xlabel '$y_2$ [m]'
set ylabel '$\dot{x_2}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 13:16 t 'Runge-Kutta 4.Ordnung', '' u 5:8 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_yvx2.eps

set output 'geschw1_yvy2.tex'
set xlabel '$y_2$ [m]'
set ylabel '$\dot{y_2}$ [m/s]'
set key outside bot center box
p 'rk2rk4_1.txt' u 13:17 t 'Runge-Kutta 4.Ordnung', '' u 5:9 t 'Runge-Kutta 2.Ordnung' 
set output
!epstopdf geschw1_yvy2.eps