
def BitMask (counter)
  return 2**counter
end

def XYZ ()
  ['X', 'Y', 'Z'].map{|direction|
    yield(direction)
  }
end

def XYZi ()
  counter = 0
  ['X', 'Y', 'Z'].map{|direction|
    ret = yield(direction, counter)
    counter+=1
    ret
  }
end

def Xyz ()
  ['x', 'y', 'z'].map{|direction|
    yield(direction, direction.upcase)
  }
end

def Xyzi ()
  counter = 0
  ['x', 'y', 'z'].map{|direction|
    ret = yield(direction, direction.upcase, counter)
    counter+=1
    ret
  }
end

def Prev(dir)
  return {'X' => 'Z', 'Y' => 'X', 'Z' => 'Y','x' => 'z', 'y' => 'x', 'z' => 'y'}[dir]
end

def Next(dir)
  return {'X' => 'Y', 'Y' => 'Z', 'Z' => 'X','x' => 'y', 'y' => 'z', 'z' => 'x'}[dir]
end

def PlusMinus(number)
  sum = [[]]
  number.times{
    tgt = []
    sum.each{|set|
      tgt << (set + ['+'])
      tgt << (set + ['-'])
    }
    sum = tgt
  }
  sum.map{|set|
    yield(set)
  }
end
