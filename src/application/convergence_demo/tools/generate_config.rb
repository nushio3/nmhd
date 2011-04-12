#!/usr/bin/env ruby

require 'nushio_script/generate_config'

class CustomConfig < Config
  def initialize(space_interpolation = :minbee2, time_integral = :midpoint2, axis = 'Z', resolution = 64, precision = 2)
    super(space_interpolation, time_integral)
    @axis = axis
    @resolution = resolution
    @precision = precision
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
    "convergence_#{precision}_#{resolution}_" + super
  end

  def to_s()
    return super + (<<END)
Normal = '#{@axis}'
Resolution = #{@resolution}
Precision = #{@precision}
END
  end
  attr_accessor :axis, :resolution, :precision
end


if __FILE__ == $0
  conf = CustomConfig::new()
  conf.parse(ARGV)
  conf.generate!()
end
