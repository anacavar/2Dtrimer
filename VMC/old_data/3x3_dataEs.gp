set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '3x3_VMC_data_variables.png'
plot 'dataEs.txt' every ::0::2 using 3:1:4 w l palette title 'gamma=4.50', 'dataEs.txt' every ::3::5 u 3:1:4 w l palette title 'gamma=4.53', 'dataEs.txt' every ::6::8 u 3:1:4 w l palette title 'gamma=4.56'
