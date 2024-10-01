unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre γ i s, uz konstantan α'
set xlabel "s / Å^{-1}"
set ylabel "E_{VMC} / K" 
set label "N_{w}=100, N_{t}=500, N_{b}=100, α=9.2 Å" at graph 0.05, 0.95 front
set terminal png size 800, 600
set output 'VMC-a-graph.png'
plot 'dataEa.txt' every ::0::9 using 5:1:4 w l palette title 'γ=4.4', 'dataEa.txt' every ::0::9 using 5:1:2 with yerrorbars notitle pointtype 0, 'dataEa.txt' every ::10::19 u 5:1:4 w l palette title 'γ=4.5', 'dataEa.txt' every ::10::19 using 5:1:2 with yerrorbars notitle pointtype 0, 'dataEa.txt' every ::20::29 u 5:1:4 w l palette title 'γ=4.6', 'dataEa.txt' every ::20::29 using 5:1:2 with yerrorbars notitle pointtype 0,'dataEa.txt' every ::30::39 using 5:1:4 w l palette title 'γ=4.7', 'dataEa.txt' every ::30::39 using 5:1:2 with yerrorbars notitle pointtype 0, 'dataEa.txt' every ::40::49 using 5:1:4 w l palette title 'γ=4.8', 'dataEa.txt' every ::40::49 using 5:1:2 with yerrorbars notitle pointtype 0
