set terminal gif animate delay 5
set output "campo_su3.gif"
stats "teste.dat" name "campo"
set xrange [campo_min_x:campo_max_x]
set yrange [-1:12]
set border 3
set tics nomirror
set grid
set xlabel "x"
set ylabel "Campo"
do for [i=0:0.1*(int(campo_blocks))] {plot "teste.dat" using 1:2 index 10*i w lines lc 7 t "", \
	"teste.dat" using 1:3 index 10*i w l lc rgb "blue" t""}
