
unset arrow
unset xrange
unset yrange
unset label
set title "DMC"
set terminal png size 1200, 600
set output 'DMC-E-graph.png'
set xlabel "blok / 500 koraka"
set ylabel "E_{DMC} / K" rotate by 90
set label "N_{w}=150, α=4.55 Å, γ=4.77, s=0.3 Å^{-1}, Δτ=0.001" at graph 0.05, 0.95 front
plot 'DMC_data.txt' u 1:2 w l title 'prosjek energije po bloku', 'DMC_data.txt' u 1:3 w l title 'prosjek energije od pocetka simulacije'

unset arrow
unset xrange
unset yrange
unset label
set title 'Promjena <r^2> kroz DMC simulaciju'
set terminal png size 1200, 600
set output 'DMC-r2-graph.png'
set xlabel "blok / 500 koraka"
set ylabel "<r^2> / Å" 
set label "N_{w}=150, α=4.55 Å, γ=4.77, s=0.3 Å^{-1}, Δτ=0.001" at graph 0.05, 0.95 front
plot 'DMC_data.txt' u 1:6 w l title 'prosjek <r^2> po bloku', 'DMC_data.txt' u 1:7 w l title 'prosjek <r^2> od pocetka simulacije'
