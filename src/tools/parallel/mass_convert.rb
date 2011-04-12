#!/usr/bin/env ruby

def sh(str)
  STDERR.puts str
  system(str)
end

ARGV.each{|dirn|
  dirs = dirn.split('/')
  body = dirs[0..-2]
  last = dirs[-1]
  outdir = (body + ['anal']).join('/')
  `mkdir -p #{outdir}`
  STDERR.puts Time.now
  sh "./tools/3.out #{dirn} #{outdir}/#{last}"
  sleep 0.1
}
