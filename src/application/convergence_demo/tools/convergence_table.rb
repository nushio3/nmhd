#!/usr/bin/env ruby

def fmt(str)
  i = str.index('.')
  return str[0..i+6]
end

(0..7).each{|logres|
  res = 10 * 2**logres
  (1..2).each{|prec|
    [1,10].each{|t|
      fn = "stat/convergence_#{prec}_#{res}_minbee2_midpoint2/dump_#{t}_energies.txt"
      mx, ekin, etot = ['-','-','-']
      if File.exist?(fn)
        str = ''
        open(fn,'r'){|fp|
          str = fp.read
        }
        words = str.split(/\s+/)
        mx = fmt(words[1])
        ekin = fmt(words[7])
        etot = fmt(words[13])
      end
      print " & #{mx} & #{ekin} & #{etot}"
    }
    print "\\\\\n"
  }
}
