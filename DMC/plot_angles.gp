set xlabel "kut u radijanima"
set ylabel "broj setaca"
set title "Distribucija kutova DMC"
set terminal png
set output 'graph_angles.png'
plot 'data_angles_DMC.txt' u 1:2 w l title 'distribucija po kutu'
