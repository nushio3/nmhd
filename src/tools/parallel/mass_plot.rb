#!/usr/bin/env ruby

Tmpfn = 'tmp.gnuplot'
def gnuplot(str)
  STDERR.puts str
  open(Tmpfn, 'w'){|fp| fp.puts str }
  `gnuplot #{Tmpfn}`
end

t=0.0
ARGV.each{|fn|
  ofn=fn[0..-5]+'.png'
  t+=0.01
  t_str = sprintf("%3.3f",t)
  gnuplot <<GNUPLOT
set term png
set out '#{ofn}'
set pm3d
set pm3d map
set size square
set xlabel 'azimuthal'
set ylabel 'radial'
set cbrange [2:3.5]
set title 'density'
splot "#{fn}" u 1:2:($4) t 't=#{t_str}sc'
GNUPLOT
}
