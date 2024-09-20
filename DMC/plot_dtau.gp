unset xrange
set title "10-6_10-5"
set terminal png
set output "10-6_10-5.png"
plot 'dataE_dtau.txt' every ::0::9 u 1:2 w l

unset xrange
set title "10-5_10-4"
set terminal png
set output "10-5_10-4.png"
plot 'dataE_dtau.txt' every ::10::19 u 1:2 w l

unset xrange
set title "10-4_10-3"
set terminal png
set output "10-4_10-3.png"
plot 'dataE_dtau.txt' every ::20::29 u 1:2 w l

unset xrange
set title "10-3_10-2"
set terminal png
set output "10-3_10-2.png"
plot 'dataE_dtau.txt' every ::30::39 u 1:2 w l

unset xrange
set title "10-2_10-1"
set terminal png
set output "10-2_10-1.png"
plot 'dataE_dtau.txt' every ::40::49 u 1:2 w l