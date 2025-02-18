set terminal png size 1920,1080 enhanced font "Times New Roman, 22"
stats "teste2.dat" name "campo"
#set xrange [campo_min_x:campo_max_x]
set xrange [campo_min_x:campo_max_x]
set yrange [campo_min_y - 0.5:10]
set border linewidth 3.0 lc rgb "black"
#set xtics linewidth 2.0
#set ytics linewidth 2.0
set tics nomirror
set grid
set xlabel "x"
set ylabel "Valor do Campo"

set output "doiscampos.png"
plot "teste2.dat" using 1:2 w l lw 3.0 lc rgb "red" title "Campo 1", \
	 "teste2.dat" using 1:3 w l lw 3.0 lc rgb "blue" title "Campo2"

