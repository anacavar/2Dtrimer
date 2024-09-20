set xlabel "r"
set ylabel ""
set terminal png
set output 'potencijal_valna.png'
set xrange[-0:15]
set yrange[-15:25]
set xzeroaxis
set yzeroaxis
set style line 1 lt 2 lc rgb "black" lw 2 dashtype 2  # Set dashed line style
set arrow from graph 0, -12 to graph 1, -12 nohead ls 1  # Dashed line at y=-12
set label "Y = -12" at graph 0.05, -12 left  # Add label to the line
plot 'output.txt' u 1:2 w l title 'Lennard-Jones potencijal', 'output.txt' u 1:3 w l title '|Ψ(x)|^2'

set xrange[-0:20]
set yrange[-0.5:0.5]
set xzeroaxis
set yzeroaxis
unset label
set output 'fdr_fddr.png'
plot 'output.txt' u 1:4 w l title 'f^{dr}(r)', 'output.txt' u 1:5 w l title 'f^{ddr}(r)'

set xrange[-0:20]
set yrange[-1:1]
set xzeroaxis
set yzeroaxis
set ylabel 'f^{dr}(r)'
set output 'fdr.png'
plot 'output.txt' u 1:4 w l

set xrange[-0:20]
set yrange[-10:10]
set xzeroaxis
set yzeroaxis
set ylabel 'f^{ddr}(r)'
set output 'fddr.png'
plot 'output.txt' u 1:5 w l
