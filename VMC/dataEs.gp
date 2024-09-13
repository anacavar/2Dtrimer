set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
set title 'VMC'
set terminal png
set output 'VMC_data_variables.png'
plot 'dataEs.txt' every ::0::4 using 3:1:4 w l palette title 'gamma=4.568', 'dataEs.txt' every ::5::9 u 3:1:4 w l palette title 'gamma=4.556', 'dataEs.txt' every ::10::14 u 3:1:4 w l palette title 'gamma=4.544', 'dataEs.txt' every ::15::19 u 3:1:4 w l palette title 'gamma=4.532', 'dataEs.txt' every ::20::24 u 3:1:4 w l palette title 'gamma=4.520'