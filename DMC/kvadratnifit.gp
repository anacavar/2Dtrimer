unset xrange
unset yrange
unset label
set title "Ovisnost DMC energije o duljini vremenskog koraka Δτ"
set xrange [0.001:0.01]
set yrange [-33.7:-33.3]
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png size 800, 600
set output "U2_dtau.png"
set label "N_{w}=100, N_{t}=500, N_{b}=100, σ=4 Å, ε=20 K" at graph 0.05, 0.95 front
plot "dataE_dtau.txt" u 1:2 w l, 'dataE_dtau.txt' using 1:2:3 with yerrorbars notitle pointtype 0


unset label
unset xrange
unset yrange
set title "Ovisnost DMC energije o duljini vremenskog koraka Δτ"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set xrange [0:0.05]
set terminal png size 800, 600
set output "U2_dtau_whole.png"
set label "N_{w}=100, N_{t}=500, N_{b}=100, σ=4 Å, ε=20 K" at graph 0.05, 0.95 front
plot "dataE_dtau.txt" u 1:2 w l

