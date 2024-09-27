
unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ"
set terminal png
set output "DMC-dtau-graph.png"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
plot 'dataE_dtau.txt' u 1:2 w l notitle

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.000001-0.000010"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-6_10-5.png"
plot 'dataE_dtau.txt' every ::0::9 u 1:2 w l notitle, 'dataE_dtau.txt' every ::0::9 using 1:2:3 with yerrorbars notitle pointtype 0

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.00001-0.00010"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-5_10-4.png"
plot 'dataE_dtau.txt' every ::10::19 u 1:2 w l notitle, 'dataE_dtau.txt' every ::0::9 using 1:2:3 with yerrorbars notitle pointtype 0

