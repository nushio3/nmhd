#!/usr/bin/env ruby

require 'nushio_script/generate_config'
require 'nushio_script/detect_machine'

class CustomConfig < Config
  def initialize(space_interpolation = :minbee2, time_integral = :midpoint2, magnetic_field = 0.1)
    config_options = {:external_force => true, :shearing_boundary => true, :external_acceleration => true}
    super(space_interpolation, time_integral, config_options)

    @parallel_mpi = true
    
    @magnetic_field = magnetic_field
    @parser.on('--magnet m', 'strength of magnetic field'){|m|
      @magnetic_field = m.to_f
    }
  end
  
  def tag ()
    "mri_B#{@magnetic_field}_" + super
  end

  def to_s()
    return super + (<<END)
# custom parameters here
MagneticField = #{@magnetic_field}
END
  end
  attr_accessor :magnetic_field
end


if __FILE__ == $0
  conf = CustomConfig::new()
  conf.parse(ARGV)
  conf.generate!()
end
