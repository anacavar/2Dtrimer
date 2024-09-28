unset xrange
unset yrange
unset label
set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'Ovisnost DMC energije o vremenskom koraku Δτ za različit broj šetača'
set ylabel "E_{DMC} / K"
set xlabel "Δτ"
set yrange [-33.9:-33.55]
set terminal png size 900, 600
set output 'DMC-dtau-Nw.png'
plot 'DATA_NW.txt' every ::0::4 using 1:2:4 w l palette title 'N_{w}=50', 'DATA_NW.txt' every ::0::4 using 1:2:3 with yerrorbars notitle pointtype 0, 'DATA_NW.txt' every ::5::9 u 1:2:4 w l palette title 'N_{w}=80', 'DATA_NW.txt' every ::5::9 using 1:2:3 with yerrorbars notitle pointtype 0, 'DATA_NW.txt' every ::10::14 u 1:2:4 w l palette title 'N_{w}=100', 'DATA_NW.txt' every ::10::14 using 1:2:3 with yerrorbars notitle pointtype 0, 'DATA_NW.txt' every ::15::19 using 1:2:4 w l palette title 'N_{w}=120', 'DATA_NW.txt' every ::15::19 using 1:2:3 with yerrorbars notitle pointtype 0, 'DATA_NW.txt' every ::20::24 using 1:2:4 w l palette title 'N_{w}=150', 'DATA_NW.txt' every ::20::24 using 1:2:3 with yerrorbars notitle pointtype 0,
