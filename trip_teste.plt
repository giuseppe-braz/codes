set terminal png size 1920,1080 enhanced font "Times New Roman, 26"
#stats "teste4.dat" name "campo"
#set xrange [campo_min_x:campo_max_x]
set output 'su2_trip.png'
set xrange [-10:10]
set yrange [0:3.2]
set border linewidth 3.0 lc rgb "black"
#set xtics linewidth 2.0
#set ytics linewidth 2.0
set tics nomirror
set grid
set xlabel "x"
set ylabel "Valor do Campo"

set key box
set key bottom right

set style line 1 lt 1 lc rgb "black" lw 2.0 pt 5 ps 1
set style line 2 lt 1 lc rgb "red" lw 2.0 pt 5 ps 1
set style line 3 lt 1 lc rgb "black" lw 2.0 pt 7 ps 1
set style line 4 lt 1 lc rgb "black" lw 2.0 pt 2 ps 1
set style line 5 lt 1 lc rgb "red" lw 2.0 pt 7 ps 1
set style line 6 lt 1 lc rgb "black" lw 2.0 pt 9 ps 1
set style line 7 lt 1 lc rgb "red" lw 2.0 pt 9 ps 1

mylist = "b=-2.0 b=-1.0 b=-0.5 b=0.0 b=0.5 b=1.0 b=2.0"

plot for [i=1:7] sprintf("trip%d.dat",i) w linespoints linestyle i t word(mylist, i)
	 
