set xlabel "r"
set ylabel ""
set terminal png
set output 'potencijal_valna.png'
set xrange[-1:15]
set yrange[-15:25]
set xzeroaxis
set yzeroaxis
plot 'output.txt' u 1:2 w l title 'Lennard-Jones potencijal', 'output.txt' u 1:3 w l title 'korelacijska funkcija'

set xrange[-1:50]
set yrange[-20:20]
set xzeroaxis
set yzeroaxis
set output 'fdr_fddr.png'
plot 'output.txt' u 1:4 w l title 'fdr', 'output.txt' u 1:5 w l title 'fddr'

