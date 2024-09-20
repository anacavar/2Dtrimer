set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre s i γ, uz konstantan α'
set xlabel "s / Å^{-1}"
set ylabel "E_{VMC} / K" rotate by 90
set label "N_{w}=300, N_{t}=1000, N_{b}=300, α=4.55" at graph 0.05, 0.95 front
set style textbox
set xrange [0.26:0.36]
set terminal png size 800,600
set output '03_BEST_7x10_VMC_data_variables_a.png'
plot '03_dataEa.txt' every ::0::9 using 5:1:4 w l palette title 'γ=4.68', '03_dataEa.txt' every ::0::9 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::10::19 u 5:1:4 w l palette title 'γ=4.72', '03_dataEa.txt' every ::10::19 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::20::29 u 5:1:4 w l palette title 'γ=4.76', '03_dataEa.txt' every ::20::29 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::30::39 u 5:1:4 w l palette title 'γ=4.80', '03_dataEa.txt' every ::30::39 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::40::49 u 5:1:4 w l palette title 'γ=4.83', '03_dataEa.txt' every ::40::49 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::50::59 u 5:1:4 w l palette title 'γ=4.87', '03_dataEa.txt' every ::50::59 using 5:1:2 with yerrorbars notitle pointtype 0, '03_dataEa.txt' every ::60::69 u 5:1:4 w l palette title 'γ=4.91', '03_dataEa.txt' every ::60::69 using 5:1:2 with yerrorbars notitle pointtype 0
