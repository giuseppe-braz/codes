set terminal png size 1920,1080 enhanced font "Times New Roman, 18"
set output "doiscampos.png"
set xrange [-15:15]
set yrange [-0.1:7]
set xlabel "x"
set ylabel "Valor do campo"
set grid
set tics nomirror

plot "teste.dat" every 5 using 1:2 with points pt 7 ps 1.3 lc rgb 'red' t 'Numérico', \
	 "teste.dat" using 1:3 w l lc 'blue' lw 3.0 t 'Analítico'

