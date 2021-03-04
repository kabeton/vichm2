set xrange [0:20]
set yrange[-1:1]
list = system('ls *.dat')
do for [file in list] {
  set terminal gif animate delay 5
  set output "gifs/".file.".gif"
  do for [i=2:36] {
    plot file u 1:i w lines  
  }
}
