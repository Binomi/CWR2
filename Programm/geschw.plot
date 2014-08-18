reset
set mulitplot
set terminal epslatex color
set output 'geschw1.tex'

set size 0.5,1.0
set origin 0.0,0.0
p x**2
set origin 0.5,0.0
p x**3

set output
unset multiplot
!epstopdf geschw1.eps