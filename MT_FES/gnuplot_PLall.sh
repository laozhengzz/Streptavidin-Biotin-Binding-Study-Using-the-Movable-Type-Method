unset key
set terminal x11
set bar 1.000000 front
set border 4095 front lt black linewidth 2.000 dashtype solid
set style circle radius graph 0.02, first 0.00000, 0.00000 
set style ellipse size graph 0.05, 0.03, first 0.00000 angle 0 units xy
set style textbox transparent margins  1.0,  1.0 border
set pm3d map
set contour base
set cntrparam levels discrete 0.0,4.0,7.0,12.0,18.0,20.0
#set cntrparam levels incremental 0.0,4.0,20.0
set key at 2,3.5,2
set dgrid3d 50 50 qnorm 3
#set palette defined (0 "royalblue",1 "turquoise",3 "yellow",5 "red")
set style line 2 lt 3 lc rgb "black" lw 10
set style line 3 lt 3 lc rgb "red" lw 10
set style line 4 lt 3 lc rgb "turquoise" lw 10
set style line 5 lt 3 lc rgb "yellow" lw 10
set style line 6 lt 3 lc rgb "green" lw 10
set style line 7 lt 3 lc rgb "orange" lw 10
set style line 8 lt 3 lc rgb "magenta" lw 10
set style line 9 lt 3 lc rgb "turquoise" lw 10
set style increment user
#set palette rgbformulae 33,13,10
#set palette rgb 3,11,6 ##good
#set palette rgb 13,1,6
set xlabel "GLY48-Ca/ILE30-Ca" font "Helvetica,24"
set xlabel offset character 4, -2, 0  
set xrange [10:20] noreverse nowriteback
set ylabel "ASN49-Ca/LEU109-CA" font "Helvetica,24"
set ylabel offset character -4, 0, 0
set yrange [14:28] noreverse nowriteback
set xtics font "Helvetica,24"
set ytics font "Helvetica,24"
set zlabel "z" 
set zlabel  offset character 10, 0, 0 font "" textcolor lt -1 norotate
set colorbox vertical origin screen 0.9,0.2,0 size screen 0.03,0.6, 0 front
splot 'grid_dPL.txt' u 1:2:3  w pm3d lw 3
#pause -1

set term postscript landscape enhanced color dashed font "Helvetica,18"
set output '|ps2pdf - PLholo.pdf'
#set lmargin at screen 0.05;
#set rmargin at screen 0.9;
#set bmargin at screen 0.1;
#set tmargin at screen 0.95;
replot
#set terminal jpeg nocrop large enhanced font "Helvetica,14" size 900,800
#set terminal jpeg nocrop large enhanced font "Helvetica,14" 
#set term postscript landscape enhanced color dashed font "Helvetica,18”
#set output "holo.ps"
#replot
