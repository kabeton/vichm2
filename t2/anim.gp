stats "anim.dat" u 1:2 nooutput
blocks = STATS_blocks

#p
set terminal gif animate delay 10
set output "gifs/p.gif"
set grid x y
show grid
set xrange[-10:10]
set xlabel "x, m"
set ylabel "P, Pa"
set xtics -10,2,10
set yrange[1e5:1e6]
do for [i=0:blocks-2] {
  plot "anim.dat" index i u 1:2 noautoscale title columnheader(1)
}

#r
set terminal gif animate delay 10
set output "gifs/r.gif"
set grid x y
show grid
set xrange[-10:10]
set xtics -10,2,10
set xlabel "x, m"
set ylabel "{/symbol r}, kg/m^3" enhanced
set yrange[0:14]
do for [i=0:blocks-2] {
  plot "anim.dat" index i u 1:3 noautoscale title columnheader(1)
}

#u
set terminal gif animate delay 10
set output "gifs/u.gif"
set grid x y
show grid
set xrange[-10:10]
set xtics -10,2,10
set xlabel "x, m"
set ylabel "u, m/s"
set yrange[0:250]
do for [i=0:blocks-2] {
  plot "anim.dat" index i u 1:4 noautoscale title columnheader(1)
}

#e
set terminal gif animate delay 10
set output "gifs/e.gif"
set grid x y
show grid
set xrange[-10:10]
set xlabel "x, m"
set ylabel "e, J/kg"
set xtics -10,2,10
set yrange[6e4:18e4]
do for [i=0:blocks-2] {
  plot "anim.dat" index i u 1:5 noautoscale title columnheader(1)
}
