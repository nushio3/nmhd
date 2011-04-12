#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'generate_config'

def run()
  sleep 0.1
  system('rake clobber')
  sleep 0.1
  system('rake')
  sleep 0.1
  system('./demo.out')
end


[[:minbee2, :midpoint2]].each{|s,t|
  ['Z'].each{|dir|
    (0..7).each{|logres|
      (1..2).each{|precision|
        res = 10*(2**logres).to_i
        STDERR.puts "resolution = #{res}; precision = #{precision};"
        CustomConfig::new(s, t, dir, res, precision).generate!()
        run()
      }
    }
  }
}

