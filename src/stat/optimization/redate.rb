#!/usr/bin/env ruby

Start = 1295096638

while line = gets
  words = line.split(/\s+/)
  t = words[0].to_i
  t -= Start
  words[0] = sprintf('%06d', t)
  puts words.join(' ')
end
