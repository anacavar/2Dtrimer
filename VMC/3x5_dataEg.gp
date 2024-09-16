set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '3x5_VMC_data_variables_g.png'
plot 'dataEs.txt' every ::0::4 using 3:1:5 w l palette title 'gamma=4.50', 'dataEs.txt' every ::5::9 u 3:1:5 w l palette title 'gamma=4.53', 'dataEs.txt' every ::10::14 u 3:1:5 w l palette title 'gamma=4.56'
