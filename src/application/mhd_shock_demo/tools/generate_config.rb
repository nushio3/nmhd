#!/usr/bin/env ruby

require 'nushio_script/generate_config'

class CustomConfig < Config
  def initialize(space_interpolation = :minbee2, time_integral = :midpoint2, axis = 'X')
    super(space_interpolation, time_integral)
    @axis = axis
    @parser.on('-axis X', '(X, Y, Z) shock normal'){|n|
      case n
      when 'X', 'Y', 'Z'
        @axis = n
      else
        throw 'unknown axis: ' + n
      end
    }
  end
  
  def tag ()
    "shock_#{@axis}_" + super
  end

  def to_s()
    return super + (<<END)
Normal = '#{@axis}'
END
  end
  attr_accessor :axis
end


if __FILE__ == $0
  conf = CustomConfig::new()
  conf.parse(ARGV)
  conf.generate!()
end
