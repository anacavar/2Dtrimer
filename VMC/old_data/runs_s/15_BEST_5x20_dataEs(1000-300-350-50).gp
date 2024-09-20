set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC za različite parametre α i γ, uz konstantan s'
set xlabel "α / Å"
set ylabel "E_{VMC} / K" rotate by 90
set label "N_{w}=300, N_{t}=1000, N_{b}=300, s=0.3" at graph 0.05, 0.95 front
set style textbox
set xrange [4.48:4.62]
set terminal png # size 800,600
set output '15_BEST_5x20_VMC_data_variables_s.png'
plot '15_dataEs.txt' every ::0::19 using 3:1:4 w l palette title 'γ=4.67', '15_dataEs.txt' every ::0::19 using 3:1:2 with yerrorbars notitle pointtype 0, '15_dataEs.txt' every ::20::39 u 3:1:4 w l palette title 'γ=4.72', '15_dataEs.txt' every ::20::39 using 3:1:2 with yerrorbars notitle pointtype 0, '15_dataEs.txt' every ::40::59 u 3:1:4 w l palette title 'γ=4.76', '15_dataEs.txt' every ::40::59 using 3:1:2 with yerrorbars notitle pointtype 0, '15_dataEs.txt' every ::60::79 using 3:1:4 w l palette title 'γ=4.81', '15_dataEs.txt' every ::60::79 using 3:1:2 with yerrorbars notitle pointtype 0, '15_dataEs.txt' every ::80::99 using 3:1:4 w l palette title 'γ=4.86','15_dataEs.txt' every ::80::99 using 3:1:2 with yerrorbars notitle pointtype 0, '15_dataEs.txt' every ::100::109 using 3:1:4 w l palette title 'γ=4.91', '15_dataEs.txt' every ::100::109 using 3:1:2 with yerrorbars notitle pointtype 0
