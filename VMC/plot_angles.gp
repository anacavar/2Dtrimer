set xlabel "kut / °"
set ylabel "vjerojatnost"
set title "Distribucija kutova VMC"
set xrange [0:180]
set yrange [0:0.004]
set arrow from 60, graph 0 to 60, graph 1 nohead lc rgb "blue"
set label "60°" at 60, 0 offset 1,1
set terminal png
set output 'VMC-angles-graph.png'
plot "data_angles.txt" u 1:2 w l title "kut α", "data_angles.txt" u 1:3 w l title "kut β", "data_angles.txt" u 1:4 w l title "kut γ"  
