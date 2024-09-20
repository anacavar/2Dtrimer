set title "Distribucija udaljenosti r_{12} i r_{13} u VMC
set pm3d map
set palette defined (0 "white", 1 "purple", 2 "blue", 3 "red", 4 "orange", 5 "yellow")
set cbrange [0:*]       # Automatically scale the color range
set xlabel "r_{12} / Å"
set ylabel "r_{13} / Å"
set xrange [0:20]
set yrange [0:20]
set colorbox vertical    # Enable color box on the right side
set key title "P(r12,r13,r23)"

# Plot the data as a heatmap
set terminal png
set output "r12_r13_graph.png"
splot 'data_r12_r13_normirano.txt' using 1:2:3 notitle with image