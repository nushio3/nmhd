#!/usr/bin/env ruby

require 'nushio_script/config_generate'

include Config

parse_config(open(ARGV[0], 'r').read).each{|d|
  puts d.key.inspect + " => " + d.value.inspect
}

