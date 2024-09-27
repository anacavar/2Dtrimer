unset xrange
unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i s, uz konstantan γ'
set xlabel "α / Å"
set ylabel "E_{VMC} / K" rotate by 90
set yrange [-34:-30]
set xrange [9:9.4]
set label "N_{w}=100, N_{t}=500, N_{b}=100, γ=7.4" at graph 0.05, 0.95 front
set terminal png
set output 'VMC-g-graph_U2.png'
plot 'dataEg.txt' every ::0::9 using 3:1:5 w l palette title 's=0.60 Å^{-1}', 'dataEg.txt' every ::0::9 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::10::19 u 3:1:5 w l palette title 's=0.62 Å^{-1}', 'dataEg.txt' every ::10::19 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::20::29 u 3:1:5 w l palette title 's=0.64 Å^{-1}', 'dataEg.txt' every ::20::29 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::30::39 using 3:1:5 w l palette title 's=0.66 Å^{-1}', 'dataEg.txt' every ::30::39 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::40::49 using 3:1:5 w l palette title 's=0.68 Å^{-1}', 'dataEg.txt' every ::40::49 using 3:1:2 with yerrorbars notitle pointtype 0
