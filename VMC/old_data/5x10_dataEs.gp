set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '5x10_VMC_data_variables.png'
plot 'dataEs.txt' every ::0::9 using 3:1:4 w l palette title 'α=4.50', 'dataEs.txt' every ::10::19 u 3:1:4 w l palette title 'α=', 'dataEs.txt' every ::20::29 u 3:1:4 w l palette title 'gamma=4.56', 'dataEs.txt' every ::30::39 using 3:1:4 w l palette title 'gamma=4.50', 'dataEs.txt' every ::40::49 using 3:1:4 w l palette title 'gamma=4.50'
