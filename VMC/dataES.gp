set palette model RGB
set palette defined ( 0 'blue', 1 'red' )
plot 'dataEs.txt' every ::0::4 using 3:1:4 w l palette, 'dataEs.txt' every ::5::9 u 3:1:4 w l palette, 'dataEs.txt' every ::10::14 u 3:1:4 w l palette, 'dataEs.txt' every ::15::19 u 3:1:4 w l palette, 'dataEs.txt' every ::20::24 u 3:1:4 w l palette