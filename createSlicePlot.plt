set term pdfcairo
set output "Parameter-Response-Slices.pdf"

set grid lc "black"

set style fill transparent solid 0.3 noborder

set log y
set format y "10^{%T}"
unset log x 

set ylabel "CPU Seconds"
set xlabel "Parameter Value

#set ytics add (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 20000, 50000)

### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
#set label 1 'a' at graph 0.92,0.9
unset yrange
unset xrange
set xlabel "keepglue"
plot "scenarios/CaDiCaL/slices-circuit-fuzz/responses/response-keepglue.csv" using 1:2 smooth unique notitle, '' using 1:2 lc 1 notitle, '' using 1:3:4 with filledcurves notitle

# --- GRAPH b
#set label 1 'b' at graph 0.92,0.9
set xlabel "BACKBONE TRIALS"
plot "scenarios/LKH/slices-tsp-rue-1000-3000/responses/response-backbone_trials.csv" using 1:2 smooth unique notitle, '' using 1:2 lc 1 notitle, '' using 1:3:4 with filledcurves notitle

# --- GRAPH c
#set label 1 'c' at graph 0.92,0.9
set xlabel "Npop"
set xrange[0:1000]
plot "scenarios/EAX/slices-tsp-rue-1000-3000/responses/response-npop.csv" using 1:2 smooth unique notitle, '' using 1:2 lc 1 notitle, '' using 1:3:4 with filledcurves notitle

# --- GRAPH d
#set label 1 'd' at graph 0.92,0.9
set xrange [120:350]
#set log x
set key bottom right
set xlabel "mip limits submipnodelim"
plot "scenarios/CPLEX/slices-regions200/responses/response-mip_limits_submipnodelim.csv" using 1:2 smooth unique notitle, '' using (0):(0) with linespoints lc 1 pt 2 title "Median of PAR10", '' using 1:2 lc 1 pt 2 notitle, '' using 1:3:4 with filledcurves lc 3 title "95% Confidence Interval"

unset multiplot
### End multiplot

