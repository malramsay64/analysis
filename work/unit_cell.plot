
set terminal png enhanced transparent truecolor size 800,800

set polar
set size ratio -1
set xrange [-10:10]
set yrange [-10:10]

unset raxis
unset rtics

set output 'unit_cell.png'

plot 'cell.dat' using 1:2:(1) with points
