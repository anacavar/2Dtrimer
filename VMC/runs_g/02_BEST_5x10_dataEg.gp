set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i s, uz konstantan γ'
set xlabel "α / Å"
set ylabel "E_{VMC} / K" rotate by 90
set xrange [4.4:4.7]
set label "N_{w}=300, N_{t}=1000, N_{b}=300, γ=4.77" at graph 0.05, 0.95 front
set terminal png
set output 'final_5x10_VMC_data_variables_g.png'
plot '02_BEST_dataEg.txt' every ::0::9 using 3:1:5 w l palette title 's=4.68', '02_BEST_dataEg.txt' every ::0::9 using 3:1:2 with yerrorbars notitle pointtype 0, '02_BEST_dataEg.txt' every ::10::19 u 3:1:5 w l palette title 's=0.29', '02_BEST_dataEg.txt' every ::10::19 using 3:1:2 with yerrorbars notitle pointtype 0, '02_BEST_dataEg.txt' every ::20::29 u 3:1:5 w l palette title 's=0.33', '02_BEST_dataEg.txt' every ::20::29 using 3:1:2 with yerrorbars notitle pointtype 0, '02_BEST_dataEg.txt' every ::30::39 using 3:1:5 w l palette title 's=0.37', '02_BEST_dataEg.txt' every ::30::39 using 3:1:2 with yerrorbars notitle pointtype 0, '02_BEST_dataEg.txt' every ::40::49 using 3:1:5 w l palette title 's=0.41', '02_BEST_dataEg.txt' every ::40::49 using 3:1:2 with yerrorbars notitle pointtype 0
