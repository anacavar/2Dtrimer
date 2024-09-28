unset xrange
unset yrange
unset label
set title "Ovisnost DMC energije o duljini vremenskog koraka Δτ"
set xrange [0.001:0.01]
set yrange [-33.7:-33.3]
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
c = -33.65
a = 4000
f(x) = a*x**2 + + c
fit f(x) "dataE_dtau.txt" u 1:2 via a, c
set terminal png size 800, 600
set output "U2-dtau-kvadratni.png"
set label "N_{w}=100, N_{t}=500, N_{b}=100, σ=4 Å, ε=20 K" at graph 0.05, 0.95 front
plot "dataE_dtau.txt" u 1:2 w l, 'dataE_dtau.txt' using 1:2:3 with yerrorbars notitle pointtype 0, f(x) title 'kvadratni fit'


unset label
unset xrange
unset yrange
set title "Ovisnost DMC energije o duljini vremenskog koraka Δτ za širok raspon veličine koraka"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set xrange [0:0.04]
set terminal png size 800, 600
set output "U2-dtau-range.png"
set label "N_{w}=100, N_{t}=500, N_{b}=100, σ=4 Å, ε=20 K" at graph 0.05, 0.95 front
plot "dataE_dtau.txt" u 1:2 w l


unset label
unset xrange
unset yrange
set title "Ovisnost DMC energije o duljini vremenskog koraka Δτ za jako male korake"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set xrange [0:0.0001]
set terminal png size 800, 600
set output "U2-dtau-detail.png"
set label "N_{w}=100, N_{t}=500, N_{b}=100, σ=4 Å, ε=20 K" at graph 0.05, 0.95 front
plot "dataE_dtau.txt" u 1:2 w l