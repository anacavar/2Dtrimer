set xlabel "kut u radijanima"
set ylabel "vjerojatnost"
set title "Distribucija kutova VMC"
set xrange [0:3.14]
set yrange [0:0.004]
set terminal png
set output 'graph_angles.png'
plot "data_angles.txt" u 1:2 w l title "kut α"
