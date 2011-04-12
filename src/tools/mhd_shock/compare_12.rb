#!/usr/bin/env ruby


fn = 'hydro'
col = 0
title = ARGV[0]
yrange = '[0:1.1]'
case title
when 'rho'
  col = 4
when 'vx'
  col = 5  
  yrange = '[-0.5:1.0]'
when 'vy'
  col = 6  
  yrange = '[-2.0:0.3]'
when 'vz'
  col = 7
  yrange = '[-0.1:0.1]'
when 'p'  
  col = 8
when 'bx'
  fn = 'bx'
  col = 4
when 'by'
   yrange = '[-1.1:1.1]'
  fn = 'by'
  col = 4
when 'bz'
  yrange = '[-0.1:0.1]'
  fn = 'bz'
  col = 4
end

puts <<GNU
set title '#{title}'
set xrange [-0.5 : 0.5]
set yrange #{yrange}
se gri
plot \
     "stat/shock_X_piecewiseconst1_euler1/dump_e_#{fn}.txt" u 1:#{col} w lp lt 2 t 'first order',\
     "stat/shock_X_minbee2_midpoint2/dump_e_#{fn}.txt" u 1:#{col} w l lt 1 lw 2 t 'second order'
GNU
