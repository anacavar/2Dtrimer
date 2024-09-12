set xlabel "udaljenost u angstremima"
set ylabel "broj setaca"
set title "DMC r12"
set terminal png
set output 'graph_r12_DMC.png'
plot 'data_r12_DMC.txt' u 1:2 w l title 'distribucija udaljenosti cestica 1 i 2'