
unset arrow
unset xrange
unset yrange
unset label
set title "VMC"
set terminal png size 1200, 600
set output 'VMC-E-graph_U2.png'
set xlabel "blok / 500 koraka"
set ylabel "E_{VMC} / K" rotate by 90
set label "N_{w}=150, α=9.2 Å, γ=7.75, s=0.75 Å^{-1}" at graph 0.05, 0.95 front
plot 'VMC_data.txt' u 1:2 w l title 'prosjek energije po bloku', 'VMC_data.txt' u 1:3 w l title 'prosjek energije od pocetka simulacije'

unset arrow
unset xrange
unset yrange
unset label
set title 'Promjena r_{12} kroz VMC simulaciju'
set terminal png size 1200, 600
set output 'VMC-r-graph_U2.png'
set xlabel "blok / 500 koraka"
set ylabel "r_{12} / Å" 
set label "N_{w}=150, α=9.2 Å, γ=7.75, s=0.75 Å^{-1}" at graph 0.05, 0.95 front
plot 'VMC_data.txt' u 1:4 w l title 'prosjek r_{12} po bloku', 'VMC_data.txt' u 1:5 w l title 'prosjek r_{12} od pocetka simulacije'

unset arrow
unset xrange
unset yrange
unset label
set title 'Promjena <r^2> kroz VMC simulaciju'
set terminal png size 1200, 600
set output 'VMC-r2-graph_U2.png'
set xlabel "blok / 500 koraka"
set ylabel "<r^2> / Å^2" 
set label "N_{w}=150, α=9.2 Å, γ=7.75, s=0.75 Å^{-1}" at graph 0.05, 0.95 front
plot 'VMC_data.txt' u 1:6 w l title 'prosjek <r^2> po bloku', 'VMC_data.txt' u 1:7 w l title 'prosjek <r^2> od pocetka simulacije'

