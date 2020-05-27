reset
set terminal pngcairo enhanced background "#ffffff" fontscale 10 size 6400, 4200 
set output 'ParetoHistLocal_95.png'

set grid
set termoption solid
set termoption enhanced

set key right bottom

set yrange [0:7000]
set xrange [0:120]
### The increment is set in seconds
set xtics
set ytics
set ylabel "Groundwater Storage Gain [TAF]"
set xlabel "Economic Cost [-]"

w = 24
ww = 28
ps = 10
#pt = 7

#set style fill solid 1 border lt -3
plot 'ParetoHistLocal_95.data' u 2:1 w points title'100' pointsize ps pointtype pt lt rgb "#a6cee3", \
     'ParetoHistLocal_95.data' u 4:3 w points title'200' pointsize ps  pointtype pt lt rgb "#1f78b4" ,\
     'ParetoHistLocal_95.data' u 6:5 w points title'300' pointsize ps  pointtype pt lt rgb "#b2df8a", \
     'ParetoHistLocal_95.data' u 8:7 w points title'400' pointsize ps  pointtype pt lt rgb "#33a02c", \
     'ParetoHistLocal_95.data' u 10:9 w points title'500' pointsize ps  pointtype pt lt rgb "#fb9a99", \
     'ParetoHistLocal_95.data' u 12:11 w points title'600' pointsize ps  pointtype pt lt rgb "#e31a1c", \
     'ParetoHistLocal_95.data' u 14:13 w points title'700' pointsize ps  pointtype pt lt rgb "#fdbf6f", \
     'ParetoHistLocal_95.data' u 16:15 w points title'800' pointsize ps  pointtype pt lt rgb "#ff7f00", \
     'ParetoHistLocal_95.data' u 18:17 w points title'900' pointsize ps  pointtype pt lt rgb "#cab2d6", \
     'ParetoHistLocal_95.data' u 20:19 w points title'1000' pointsize ps  pointtype pt lt rgb "#6a3d9a"



### Run from terminal
###  load 'Local_90.plt'

