reset
set terminal pngcairo enhanced background "#ffffff" fontscale 10 size 6400, 4200 
set output 'HyperVolumeLocal_95.png'

set grid
set termoption solid
set termoption enhanced

set key right bottom

set yrange [1:17]
set y2range [0:60]
set xrange [0:1000]
### The increment is set in seconds
set xtics 100, 100,1000
set ytics
set xlabel "Generation"
set ylabel "Hypervolume [-]"
set y2label "# of pareto solutions"

w = 32
ww = 28

#set style fill solid 1 border lt -3
plot 'HyperVolumeLocal_95.data' u 1:2 w line title "Hypervolume" lw w lt rgb "#1f78b4" axes x1y1, \
     'HyperVolumeLocal_95.data' u 1:3 w line title "# solutions" lw w lt rgb "#e31a1c" axes x1y2, \


### Run from terminal
###  load 'Local_90.plt'

