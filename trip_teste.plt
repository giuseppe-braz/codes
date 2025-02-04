set terminal png size 1920,1080 enhanced font "Times New Roman, 22"
#stats "teste4.dat" name "campo"
#set xrange [campo_min_x:campo_max_x]
set output 'trip_teste.png'
set xrange [-10:10]
set yrange [0:3.5]
set border linewidth 3.0 lc rgb "black"
#set xtics linewidth 2.0
#set ytics linewidth 2.0
set tics nomirror
set grid
set xlabel "x"
set ylabel "Valor do Campo"

set style line 1 lc rgb 'black' lt 1 lw 2 pt 0
set style line 2 lc rgb 'black' lt 2 lw 2 pt 0
set style line 3 lc rgb 'black' lt 1 lw 2 pt 7
set style line 4 lc rgb 'black' lt 2 lw 2 pt 7
set style line 5 lc rgb 'black' lt 1 lw 2 pt 13
set style line 6 lc rgb 'black' lt 2 lw 2 pt 13
set style line 7 lc rgb 'black' lt 1 lw 2 pt 15





plot for [i=1:7] sprintf("trip%d.dat",i) w l ls i  t ''
	 
