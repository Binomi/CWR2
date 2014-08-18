reset
set terminal epslatex color
set output 'trajektorie1.tex'
set xlabel '$x$ [m]'
set ylabel '$y$ [m]'
set key outside bot center box

set size ratio -1

p 'rk2rk4_1.txt' u 12:13 w lp t 'Runge-Kutta 4.Ordnung', '' u 4:5 t 'Runge-Kutta 2.Ordnung' w lp

set output
!epstopdf trajektorie1.eps

set output 'trajektorie2.tex'
p 'rk2rk4_2.txt' u 12:13 w lp t 'Runge-Kutta 4.Ordnung', '' u 4:5 t 'Runge-Kutta 2.Ordnung' w lp
set output
!epstopdf trajektorie2.eps

set output 'trajektorie3.tex'
p 'rk2rk4_3.txt' u 12:13 w lp t 'Runge-Kutta 4.Ordnung', '' u 4:5 t 'Runge-Kutta 2.Ordnung' w lp
set output
!epstopdf trajektorie3.eps

set output 'trajektorie4.tex'
p 'rk2rk4_4.txt' u 12:13 w lp t 'Runge-Kutta 4.Ordnung', '' u 4:5 t 'Runge-Kutta 2.Ordnung' w lp
set output
!epstopdf trajektorie4.eps

set output 'trajektorie5.tex'
p 'rk2rk4_5.txt' u 12:13 w lp t 'Runge-Kutta 4.Ordnung', '' u 4:5 t 'Runge-Kutta 2.Ordnung' w lp
set output
!epstopdf trajektorie5.eps