#!/usr/bin/gnuplot

set title ""
set logscale
set terminal pngcairo size 350,265 enhanced font 'Verdana,10'
set output '../../results/25/mae_per_atom.png'
set format y "10^%T"
set format x "10^%T"
set style line 1 lc rgb 'red' lt 1 lw 1.5 pt 1 ps 1
set style line 2 lc rgb 'green' lt 1 lw 1.5 pt 2 ps 1
set style line 3 lc rgb 'blue' lt 1 lw 1.5 pt 3 ps 1
set style line 4 lc rgb 'orange' lt 1 lw 1.5 pt 4 ps 1
set style line 5 lc rgb '#b01e7a' lt 1 lw 1.5 pt 5 ps 1
set border linewidth 1.5
#set key at 50,112
#set style line 1 lc rgb '#0060ad' lt 1 lw 1.5 pt 7 pi -1 ps 1
#set style line 2 lc rgb '#dd181f' lt 1 lw 1.5 pt 7 pi -1 ps 1
set ylabel "MAE Per Atom"
set xlabel "Training Set Size"

plot '../../results/25/phase_1/mae.dat' u 1:($3*1e3) title 'Phase 1' with linespoints ls 1, '../../results/25/phase_2/mae.dat' u 1:($3*1e3) title 'Phase 2' with linespoints ls 2, '../../results/25/phase_3/mae.dat' u 1:($3*1e3)  title 'Phase 3' with linespoints ls 3, '../../results/25/phase_all/mae.dat' u 1:($3*1e3)  title 'Phase All' with linespoints ls 4, '../../results/25/felix/mae.dat' u 1:($3*1e3)  title 'IJQC' with linespoints ls 5

