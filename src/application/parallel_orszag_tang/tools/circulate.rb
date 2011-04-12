#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'generate_config'

def sh(str)
  STDERR.puts str
  system(str)
end

[[:minbee2, :midpoint2], [:piecewiseconst1, :euler1]].each{|s,t|
  ['Y', 'Z', 'X'].each{|dir|
    sleep(0.1)
    STDERR.puts "**** testcase #{s} #{t} #{dir} ****"
    
    sh('rake clobber')
    CustomConfig::new(s, t, dir).generate!()
    sh('rake')
    sh('mpirun -n 1 ./demo.out --new')
    #sh('mpirun -n 3 ./demo.out --continue')
    #sh('mpirun -n 3 ./demo.out --continue')
  }
}

