set output "animation_free_particle.gif"
set terminal gif animate delay 20
stats "data_free_particle.dat" name "A"
set xrange [A_min_x:A_max_x]
set yrange [A_min_y-1:A_max_y+1]
set grid
set xlabel "x"
set ylabel "probability density"
do for [i=0:int(A_blocks-1)] {plot "data_free_particle.dat" index i w lines title "time ".(i+1)}
