#!/usr/bin/env ruby

def gnuplot(str)
  open('tmp','w'){|fp| fp.puts str}
  `gnuplot tmp`
end

dirs = `ls -d1 stat-mri/mri_B*`.split(/\n/).map{|s|s.strip}


dirs.each{|dir|
  title = dir.split('/')[1].split('_')[1]
  cmd = "\"< grep 'By_sq' #{dir}/log/0000.txt | awk '{print $2, sqrt($6)/$4*sqrt(2)*100}'\" t '#{title}'"
  logcmd = "\"< grep 'By_sq' #{dir}/log/0000.txt | awk '{print log($2), log(sqrt($6)/$4*sqrt(2)*100)}'\" t '#{title}'"
  
  ext = '.png'
  term = 'set term png'
  
  gnuplot <<GNU
#{term}
set xlabel 'time'
set ylabel 'perturbation amplitude'
set out 'jitujitu_#{title}#{ext}'
set yrange [0:15]
set grid
plot #{cmd}
GNU

  gnuplot <<GNU
#{term}
set xlabel 'time'
set ylabel 'perturbation amplitude'
set out 'loglog_#{title}#{ext}'
set log xy
set grid
plot #{cmd}
GNU

  gnuplot <<GNU
#{term}
set xlabel 'log(time)'
set ylabel 'log(perturbation amplitude)'
set out 'loglogzoom_#{title}#{ext}'
set xrange [0:2]
set grid
plot #{logcmd}
GNU


}
