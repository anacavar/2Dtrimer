unset arrow
set title 'VMC'
set terminal png size 1200, 600
set output 'VMC.png'
set xlabel "blok / 1000 koraka"
set ylabel "E_{VMC} / K" rotate by 90
set arrow 1 from 0,-5.34 to 1000,-5.34 nohead lc rgb "#D3D3D3" lw 1 dt 1
set label "N_{w}=500, N_{t}=1000, N_{b}=1000, α=4.55, γ=4.77, s=0.30" at graph 0.05, 0.95 front
plot 'data.txt' u 1:2 w l title 'prosjek energije po bloku', 'data.txt' u 1:3 w l title 'prosjek energije od pocetka simulacije'