set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '3x5_VMC_data_variables_g.png'
plot 'dataEg.txt' every ::0::4 using 3:1:5 w l palette title 's=0.20', 'dataEg.txt' every ::5::9 u 3:1:5 w l palette title 's=0.27', 'dataEg.txt' every ::10::14 u 3:1:5 w l palette title 's=0.33'
