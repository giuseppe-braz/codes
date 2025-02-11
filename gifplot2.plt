set terminal gif animate delay 5
set output "gif_campo.gif"
stats "teste4.dat" name "campo"
set xrange [campo_min_x:campo_max_x]
set yrange [-1:13]
set border 3
set tics nomirror
set grid
set xlabel "x"
set ylabel "Campo"
do for [i=0:(int(campo_blocks)*0.1)] {plot "teste4.dat" index 10*i w lines lc 7 t ""}
