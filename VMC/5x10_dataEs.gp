unset xrange
unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i γ, uz konstantan s'
set xlabel "γ"
set ylabel "E_{VMC} / K"
set label "N_{w}=100, N_{t}=500, N_{b}=100, s=0.68 Å^{-1}" at graph 0.05, 0.95 front
set style textbox
set yrange [-5:0]
set terminal png size 800, 600
set output 'VMC-s-graph.png'
plot 'dataEs.txt' every ::0::9 using 4:1:3 w l palette title 'α=4.1 Å', 'dataEs.txt' every ::0::9 using 4:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::10::19 u 4:1:3 w l palette title 'α=4.3 Å', 'dataEs.txt' every ::10::19 using 4:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::20::29 u 4:1:3 w l palette title 'α=4.5 Å', 'dataEs.txt' every ::20::29 using 4:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::30::39 using 4:1:3 w l palette title 'α=4.7 Å', 'dataEs.txt' every ::30::39 using 4:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::40::49 using 4:1:3 w l palette title 'α=4.9 Å','dataEs.txt' every ::40::49 using 4:1:2 with yerrorbars notitle pointtype 0
