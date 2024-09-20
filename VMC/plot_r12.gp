set xlabel "r_{ij} / Å"
set ylabel "vjerojatnost"
set title "Distribucija udaljenosti r_{ij} u VMC"
set yrange [0:0.0065]
set terminal png
set output 'graph_r12.png'
set arrow from 6.5091, graph 0 to 6.5091, graph 1 nohead lc rgb "blue"
set label "r_{RMS} = 6.5091" at 6.5091, 0 offset 1,1
plot 'data_r12.txt' u 1:2 w l title 'distribucija udaljenosti cestica 1 i 2'