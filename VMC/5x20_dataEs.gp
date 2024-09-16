set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output '5x20_VMC_data_variables.png'
plot 'dataEs.txt' every ::0::19 using 3:1:4 w l palette title 'gamma=4.50', 'dataEs.txt' every ::20::39 u 3:1:4 w l palette title 'gamma=4.53', 'dataEs.txt' every ::40::59 u 3:1:4 w l palette title 'gamma=4.56', 'dataEs.txt' every ::60::79 using 3:1:4 w l palette title 'gamma=4.50', 'dataEs.txt' every ::80::99 using 3:1:4 w l palette title 'gamma=4.50'
