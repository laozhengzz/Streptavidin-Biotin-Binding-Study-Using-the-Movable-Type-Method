module load MATLAB-compiler/R2016a
./plot_grid 0.5
gnuplot gnuplot_Pall.sh 
gnuplot gnuplot_PLall.sh 

rm grid_dP.pmf;
rm grid_dPL.pmf;
