set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '5x10_VMC_data_variables_g.png'
plot 'dataEg.txt' every ::0::9 using 3:1:5 w l palette title 's=0.25', 'dataEg.txt' every ::10::19 u 3:1:5 w l palette title 's=0.29', 'dataEg.txt' every ::20::29 u 3:1:5 w l palette title 's=0.33', 'dataEg.txt' every ::30::39 using 3:1:5 w l palette title 's=0.37', 'dataEg.txt' every ::40::49 using 3:1:5 w l palette title 's=0.41'
