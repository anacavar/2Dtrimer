unset xrange
unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i γ, uz konstantan s'
set xlabel "α / Å"
set ylabel "E_{VMC} / K"
set label "N_{w}=100, N_{t}=1000, N_{b}=200, s=0.2 Å^{-1}" at graph 0.05, 0.95 front
set style textbox
set terminal png 
set output 'VMC-s-graph_U2.png'
plot 'dataEs.txt' every ::0::4 using 3:1:4 w l palette title 'γ=11.00', 'dataEs.txt' every ::0::4 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::5::9 u 3:1:4 w l palette title 'γ=12.33', 'dataEs.txt' every ::5::9 using 3:1:2 with yerrorbars notitle pointtype 0, 'dataEs.txt' every ::10::14 u 3:1:4 w l palette title 'γ=13.66', 'dataEs.txt' every ::10::14 using 3:1:2 with yerrorbars notitle pointtype 0
