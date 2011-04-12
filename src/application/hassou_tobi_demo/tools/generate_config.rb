#!/usr/bin/env ruby

require 'nushio_script/generate_config'
require 'nushio_script/detect_machine'

class CustomConfig < Config
  def initialize(space_interpolation = :minbee2, time_integral = :midpoint2, axis = 'Z')
    config_options = {:external_force => true, :shearing_boundary => true}
    super(space_interpolation, time_integral, config_options)

    @parallel_mpi = true
    
    @axis = axis
    @parser.on('--axis X', '(X, Y, Z) shock normal'){|n|
      case n
      when 'X', 'Y', 'Z'
        @axis = n
      else
        throw 'unknown axis: ' + n
      end
    }
  end
  
  def tag ()
    "orszag_tang_#{@axis}_" + super
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
