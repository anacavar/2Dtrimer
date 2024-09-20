set xlabel "udaljenost u angstremima"
set ylabel "broj setaca"
set title "VMC r12"
set label "N_{w}=300, N_{t}=1000, N_{b}=300, α=4.55, γ=4.77, s=0.30" at graph 0.05, 0.95 front
set terminal png
set output 'graph_r12.png'
plot 'data_r12.txt' u 1:2 w l title 'distribucija udaljenosti cestica 1 i 2'