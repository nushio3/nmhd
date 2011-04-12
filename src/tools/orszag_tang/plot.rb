#!/usr/bin/env ruby
puts <<GNUPLOT
set xrange [0:#{2*Math::PI}]
set yrange [0:#{2*Math::PI}]
set pm3d map
set pm3d
se gri
splot "stat/orszag_tang_Z_minbee2_midpoint2/dump_e_hydro.txt" u 1:2:($8/$4)
GNUPLOT
