set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '5x10_VMC_data_variables_a.png'
plot 'dataEa.txt' every ::0::9 using 5:1:4 w l palette title 'gamma=4.50', 'dataEa.txt' every ::10::19 u 5:1:4 w l palette title 'gamma=4.53', 'dataEa.txt' every ::20::29 u 5:1:4 w l palette title 'gamma=4.56', 'dataEa.txt' every ::30::39 using 5:1:4 w l palette title 'gamma=4.50', 'dataEa.txt' every ::40::49 using 5:1:4 w l palette title 'gamma=4.50'
