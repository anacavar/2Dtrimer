unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre γ i s, uz konstantan α'
set xlabel "s / Å^{-1}"
set ylabel "E_{VMC} / K" 
set label "N_{w}=100, N_{t}=1000, N_{b}=200, α=8.3 Å" at graph 0.05, 0.95 front
set terminal png
set output 'VMC-a-graph_U2.png'
plot 'dataEa.txt' every ::0::4 using 5:1:4 w l palette title 'γ=12.00', 'dataEa.txt' every ::0::4 using 5:1:2 with yerrorbars notitle pointtype 0, 'dataEa.txt' every ::5::9 u 5:1:4 w l palette title 'γ=13.17', 'dataEa.txt' every ::5::9 using 5:1:2 with yerrorbars notitle pointtype 0, 'dataEa.txt' every ::10::14 u 5:1:4 w l palette title 'γ=14.33', 'dataEa.txt' every ::10::14 using 5:1:2 with yerrorbars notitle pointtype 0
