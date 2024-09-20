
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
plot 'dataE_dtau.txt' every ::0::9 u 1:2 w l notitle

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.00001-0.00010"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-5_10-4.png"
plot 'dataE_dtau.txt' every ::10::19 u 1:2 w l notitle

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.0001-0.0010"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-4_10-3.png"
plot 'dataE_dtau.txt' every ::20::29 u 1:2 w l notitle

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.001-0.010"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-3_10-2.png"
plot 'dataE_dtau.txt' every ::30::39 u 1:2 w l notitle

unset xrange
set title "Ovisnost energije DMC o vremenskom koraku Δτ na rasponu 0.01-0.10"
set xlabel "Δτ"
set ylabel "E_{DMC} / K"
set terminal png
set output "10-2_10-1.png"
plot 'dataE_dtau.txt' every ::40::49 u 1:2 w l notitle