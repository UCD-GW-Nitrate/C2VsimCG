reset
set terminal pngcairo enhanced background "#ffffff" fontscale 8 size 6400, 3200 
set output 'DTS_local_90.png'

set grid
set termoption solid
set termoption enhanced
set xdata time
set timefmt "%Y-%m-%d"
set format x "%Y"

# Legend
set key left top

set title "Excess flow above 90^{th} percentile"
show title

set yrange [0:7]
set xrange ["1965-09-01":"2009-09-01"]
### The increment is set in seconds
set xtics "1920-10-01", 315360000,"2010-10-01"

set xlabel "Time"
set ylabel "MAF"

w = 18

#set style fill solid 1 border lt -3
plot 'dts_local_90.data' u 1:2 w line title 'Friant-Kern' lw w lt rgb "#377eb8", \
     'dts_local_90.data' u 1:3 w line title 'Kern' lw w lt rgb "#4daf4a", \
     'dts_local_90.data' u 1:4 w line title 'Kings' lw w lt rgb "#984ea3", \
     'dts_local_90.data' u 1:5 w line title 'Kaweah' lw w lt rgb "#ff7f00", \
     'dts_local_90.data' u 1:6 w line title 'Tule' lw w lt rgb "#cccc29"

### Run from terminal
###  load 'Local_90.plt'

