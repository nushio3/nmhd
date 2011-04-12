#!/usr/bin/env ruby

$:.unshift(File.dirname(__FILE__))

require 'generate_config'

def run()
  system('rake clobber')
  system('rake')
  system('./demo.out')
end


[[:piecewiseconst1, :euler1], [:minbee2, :midpoint2]].each{|s,t|
  ['Y', 'Z', 'X'].each{|dir|
    CustomConfig::new(s, t, dir).generate!()
    run()
  }
}

