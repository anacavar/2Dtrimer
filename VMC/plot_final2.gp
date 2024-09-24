unset label
unset key
unset size
unset xrange
unset yrange
set xlabel "X_{E}"
set ylabel "Y_{ρ}" 
set label "σ=4Å, m=4u" at graph 0.05, 0.95 front
set key
set logscale x 10
set logscale y 10
plot "final_data.txt" using 1:2 w l notitle, "final_data.txt" every ::0::0 using 1:2 with points pt 7 title "ε = 4K" ,"final_data.txt" every ::1::1 using 1:2 with points pt 7 title "ε = 6K", "final_data.txt" every ::2::2 using 1:2 with points pt 7 title "ε = 8K", "final_data.txt" every ::3::3 using 1:2 with points pt 7 title "ε = 12K", "final_data.txt" every ::4::4 using 1:2 with points pt 7 title "ε = 16.4K", "final_data.txt" every ::5::5 using 1:2 with points pt 7 title "ε = 24.8K", "final_data.txt" every ::6::6 using 1:2 with points pt 7 title "ε = 33.2K", "final_data.txt" every ::7::7 using 1:2 with points pt 7 title "ε = 41.6K", "final_data.txt" every ::8::8 using 1:2 with points pt 7 title "ε = 50K", 1/(4*x)**(1.0/3) w l title 'klasična linija'

