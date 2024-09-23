unset xrange
unset yrange
unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i s, uz konstantan γ'
set xlabel "α / Å"
set ylabel "E_{VMC} / K" rotate by 90
set label "N_{w}=100, N_{t}=1000, N_{b}=200, γ=12" at graph 0.05, 0.95 front
set terminal png
set output 'VMC-g-graph_U2.png'
plot 'dataEg.txt' every ::0::4 using 3:1:5 w l palette title 's=0.35', 'dataEg.txt' every ::0::4 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::5::9 u 3:1:5 w l palette title 's=0.42', 'dataEg.txt' every ::5::9 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEg.txt' every ::10::14 u 3:1:5 w l palette title 's=0.48', 'dataEg.txt' every ::10::14 using 3:1:2 with yerrorbars notitle pointtype 0
