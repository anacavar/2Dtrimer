set xlabel "udaljenost u angstremima"
set ylabel "broj setaca"
set title "VMC r12"
set terminal png
set output 'graph_r12.png'
plot 'data_r12.txt' u 1:2 w l title 'distribucija udaljenosti cestica 1 i 2'