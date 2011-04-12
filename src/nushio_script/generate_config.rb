#!/usr/bin/env ruby

require 'optparse'

CONFIG_FN = 'simulation_config_extern.rb'

class Config
  attr_accessor :ready,  :space_interpolation, :time_integral
  attr_accessor :parallel_mpi, :config_options
  attr_accessor :parser

  def initialize(space_interpolation = :minbee2, time_integral = :midpoint2, config_options = {})
    @ready = true
    @space_interpolation = space_interpolation
    @time_integral = time_integral
    @config_options  = config_options 
    @parallel_mpi = false
    @parser = OptionParser.new

    @parser.banner = "generate #{CONFIG_FN}"
    @parser.on('-s', '--space=N',  '(1 or 2) the order of space interpolation'){|n|
      case n.to_i
      when 1
        @space_interpolation = :piecewiseconst1
      when 2
        @space_interpolation = :minbee2
      end
    }

    @parser.on('-t', '--time=N', '(1 or 2) the order of time integral'){|n|
      case n.to_i
      when 1
        @time_integral = :euler1
      when 2
        @time_integral = :midpoint2
      end
    }

    @parser.on('--default', 'generate config file instead of showing help when no argument is given') {|n|
    }
  end

  def tag ()
    "#{space_interpolation}_#{time_integral}"
  end
  
  def to_s()
    return <<END
SpaceInterpolation = :#{space_interpolation}
TimeIntegral = :#{time_integral}
SimulationTag = '#{tag()}'

UseMpi = #{@parallel_mpi}

ConfigOptions = #{config_options.inspect}

END
  end

  def parse(argv)
    @ready = (argv.length > 0)
    return @parser.parse(argv)
  end
  def generate!()
    unless @ready
      STDERR.puts @parser
      return false
    end
    open(CONFIG_FN,'w') {|fp|
      fp.puts to_s()
    }
    return true
  end
  
end


  
