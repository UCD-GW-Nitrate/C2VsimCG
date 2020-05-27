reset
set terminal pngcairo enhanced background "#ffffff" fontscale 6 size 6400, 3200 
set output 'FriantKernHyd.png'

set grid
set termoption solid
set termoption enhanced
set xdata time
set timefmt "%Y-%m-%d"
set format x "%Y"

# Legend
set key center top

set title "Caclulation of excess flow for Friant-Kern diversion"
show title

set yrange [0:0.9]
set xrange ["1921-10-01":"2009-09-01"]
set y2range [0:16]
### The increment is set in seconds
set xtics "1920-10-01", 315360000,"2010-10-01"
set ytics
set y2tics
set xlabel "Time"
set ylabel "Monthly water flow [MAF]"
set y2label "Cumulative water flow [MAF]"

w = 18
ww = 28
p90(x) = 0.2039649
p95(x) = 0.3258237

#set style fill solid 1 border lt -3
plot 'FrianKernHyd.data' u 1:2 w line title 'Hydrograph' lw w lt rgb "#377eb8" axes x1y1, \
     p90(x) w line title "90^{th} perc" lw 16 lt rgb "#e41a1c" dt 2, \
     p95(x) w line title "95^{th} perc" lw 16 lt rgb "#4daf4a" dt 2, \
     'FrianKernHyd.data' u 1:3 w line title 'Cumulative flow above 90^{th}%' lw ww lt rgb "#e41a1c" axes x1y2, \
     'FrianKernHyd.data' u 1:4 w line title 'Cumulative flow above 95^{th}%' lw ww lt rgb "#4daf4a" axes x1y2


### Run from terminal
###  load 'Local_90.plt'

