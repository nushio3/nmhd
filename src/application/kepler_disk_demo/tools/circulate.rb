#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'generate_config'

def sh(str)
  STDERR.puts str
  system(str)
end

magnets = (0..20).to_a.map{|p| 10**(-2+p/10.0) }
puts magnets.inspect


[[:minbee2, :midpoint2]].each{|s,t|
  magnets.each{|b_z|
    sleep(0.1)
    sh('rake clobber')
    CustomConfig::new(s, t, b_z).generate!()
    sh('rake')
    STDERR.puts "**** testcase #{s} #{t} #{b_z} ****"
    sh('mpirun -n 3 ./demo.out --new')
  }
}

