set  autoscale                        # scale axes automatically
set datafile separator ","
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
# set style fill transparent solid 0.5 noborder
set title "Cluster 2D Data Plot"
set pm3d map
set key off

splot "gkfcm-gnuplot.dat" using 1:2:3 with points lt pal pt 6
pause -1 "\nPush 'q' and 'return' to exit Gnuplot ...\n"
