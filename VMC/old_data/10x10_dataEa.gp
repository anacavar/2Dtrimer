set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '10x10_VMC_data_variables_a.png'
plot 'dataEa.txt' every ::0::9 using 5:1:4 w l palette title 'gamma=4.68', 'dataEa.txt' every ::10::19 u 5:1:4 w l palette title 'gamma=4.72', 'dataEa.txt' every ::20::29 u 5:1:4 w l palette title 'gamma=4.76', 'dataEa.txt' every ::30::39 using 5:1:4 w l palette title 'gamma=4.80', 'dataEa.txt' every ::40::49 using 5:1:4 w l palette title 'gamma=4.83', 'dataEa.txt' every ::50::59 using 5:1:4 w l palette title 'gamma=4.87', 'dataEa.txt' every ::60::69 using 5:1:4 w l palette title 'gamma=4.91', 'dataEa.txt' every ::70::79 using 5:1:4 w l palette title 'gamma=4.91', 'dataEa.txt' every ::80::89 using 5:1:4 w l palette title 'gamma=4.91', 'dataEa.txt' every ::90::99 using 5:1:4 w l palette title 'gamma=4.91'
