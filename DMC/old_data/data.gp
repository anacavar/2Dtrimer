
unset arrow
unset xrange
unset yrange
unset label
set title "DMC"
set terminal png size 1200, 600
set output 'DMC-E-graph.png'
set xlabel "blok / 1000 koraka"
set ylabel "E_{DMC} / K" rotate by 90
set label "N_{w}=500, N_{t}=1000, N_{b}=500, α=4.55, γ=4.77, s=0.30, Δτ=0.01" at graph 0.05, 0.95 front
plot 'data.txt' u 1:2 w l title 'prosjek energije po bloku', 'data.txt' u 1:3 w l title 'prosjek energije od pocetka simulacije'

unset arrow
unset xrange
unset yrange
unset label
set title 'Promjena <r^2> kroz DMC simulaciju'
set terminal png size 1200, 600
set output 'DMC-r2-graph.png'
set xlabel "blok / 1000 koraka"
set ylabel "<r^2> / Å" 
set label "N_{w}=500, N_{t}=1000, N_{b}=500, α=4.55, γ=4.77, s=0.30" at graph 0.05, 0.95 front
plot 'data.txt' u 1:6 w l title 'prosjek <r^2> po bloku', 'data.txt' u 1:7 w l title 'prosjek <r^2> od pocetka simulacije'
