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
set terminal png size 800, 600
set output "VMC-X-Y.png"
plot "data_rescaled.txt" using 1:2 w l notitle, "data_rescaled.txt" every ::0::0 using 1:2 with points pt 7 title "ε = 4K" ,"data_rescaled.txt" every ::1::1 using 1:2 with points pt 7 title "ε = 6K", "data_rescaled.txt" every ::2::2 using 1:2 with points pt 7 title "ε = 8K", "data_rescaled.txt" every ::3::3 using 1:2 with points pt 7 title "ε = 12K", "data_rescaled.txt" every ::4::4 using 1:2 with points pt 7 title "ε = 16.4K", "data_rescaled.txt" every ::5::5 using 1:2 with points pt 7 title "ε = 24.8K", "data_rescaled.txt" every ::6::6 using 1:2 with points pt 7 title "ε = 33.2K", "data_rescaled.txt" every ::7::7 using 1:2 with points pt 7 title "ε = 41.6K", "data_rescaled.txt" every ::8::8 using 1:2 with points pt 7 title "ε = 50K", 1.02*(3./(2*x))**(1.0/3) w l title 'klasična linija'

