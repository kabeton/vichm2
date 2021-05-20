stats "p.dat" u 1:2 nooutput
blocks = STATS_blocks

#p
set terminal gif animate delay 2 size 1300, 800
set output "p.gif"
set grid x y
show grid
set xrange[0:500]
set xlabel "x, m"
set ylabel "P, Pa"
set yrange[0.5e7:1.5e7]
do for [i=0:blocks-2] {
  plot "p.dat" index i u 1:2 w lp lw 2 pt 7 noautoscale title columnheader(1)
}