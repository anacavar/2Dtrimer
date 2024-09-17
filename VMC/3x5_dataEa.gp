set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '3x5_VMC_data_variables_a.png'
plot 'dataEa.txt' every ::0::4 using 5:1:4 w l palette title 's=0.20', 'dataEa.txt' every ::5::9 u 5:1:4 w l palette title 's=0.27', 'dataEa.txt' every ::10::14 u 5:1:4 w l palette title 's=0.33'
