unset label
unset arrow
set xlabel "r_{12} / Å"
set ylabel "vjerojatnost"
set title "Distribucija udaljenosti r_{12} u DMC"
set xrange [0:20]
set yrange [0:0.008]
set terminal png
set output 'DMC-r12-graph.png'
set arrow from 6.72, graph 0 to 6.72, graph 1 nohead lc rgb "blue"
set label "r_{RMS} = 6.72" at 6.72, 0 offset 1,1
set xrange [0:20]
plot 'data_r12_DMC.txt' u 1:2 w l title 'distribucija udaljenosti cestica 1 i 2'