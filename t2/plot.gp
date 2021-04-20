set terminal png
set grid x y
show grid

set output "pics/p15.png"
plot "t15.dat" u 1:2 

set output "pics/rho15.png"
plot "t15.dat" u 1:3 

set output "pics/u15.png"
plot "t15.dat" u 1:4 

set output "pics/e15.png"
plot "t15.dat" u 1:5 

set output "pics/p2.png"
plot "t2.dat" u 1:2 

set output "pics/rho2.png"
plot "t2.dat" u 1:3 

set output "pics/u2.png"
plot "t2.dat" u 1:4 

set output "pics/e2.png"
plot "t2.dat" u 1:5 