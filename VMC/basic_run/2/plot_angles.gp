set xlabel "kut u radijanima"
set ylabel "broj setaca"
set title "Distribucija kutova VMC"
set terminal png
set output 'graph_angles.png'
plot 'data_angles.txt' u 1:2 w l title 'distribucija po kutu'
