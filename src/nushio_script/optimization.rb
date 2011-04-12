module Optimization
  CanonFilename = 'canonical_result.bin'

  Unknown = 0

  class EvolutionInstruction
    def initialize(bro = 0, dra = 0.5)
      @broadness = bro
      @drasticness = [[dra, 0.0].max, 1.0].min
    end
    attr_accessor :broadness, :drasticness
  end
  
  class OptimizationSetting
    def initialize(inst = EvolutionInstruction::new())
      @choice_hash = {}
      @tag_to_be_evolved = []
      @evolution_instruction = inst
      @score = Unknown
    end
    
    def merge(other_setting)
      other = other_setting.choice_hash
      ks = choice_hash.keys.clone
      ks.each{|k|
        choice_hash[k] = other[k] if rand() < 0.5 && other.key?(k)
      }
    end

    attr_accessor :choice_hash, :tag_to_be_evolved, :evolution_instruction, :score
  end

  def self.optimize() 
    opt = Optimization::new()
    yield(opt)
    opt.write()
  end
  
  InitialFilename = 'opt_initial.bin'
  CurrentFilename = 'opt_current.bin'

  class Optimization
    def initialize()
      unless File.exist?(InitialFilename)
        @initial_setting = OptimizationSetting::new()
        STDERR.puts "EVOINST none"
      else
        open(InitialFilename, 'r'){|fp|
          @initial_setting = Marshal.load(fp)
        }
        STDERR.puts "EVOINST #{@initial_setting.tag_to_be_evolved.inspect}"
      end

      unless File.exist?(CurrentFilename)
        @current_setting = @initial_setting
      else
        open(CurrentFilename, 'r'){|fp|
          @current_setting = Marshal.load(fp)
        }
      end
    end
    
    # don't forget to call this at the end of your program, or use the optimize() block.
    def write()
      open(CurrentFilename, 'w'){|fp|
        Marshal.dump(@current_setting, fp)
      }
    end

    # most general evolution choice.
    # 1. give a string tag to distinguish choices.
    # 2. give a block that, given an evolution instruction and initial value,
    #                       returns a new value.
    def choice(tag, &evolver)
      new_value = nil
      if (not @initial_setting) || (initial_value = @initial_setting.choice_hash[tag]).nil?
        new_value = evolver.call(nil, nil)
      elsif @initial_setting.tag_to_be_evolved.index(tag)
        new_value = evolver.call(@initial_setting.evolution_instruction, initial_value)
      else
        new_value = initial_value
      end
      @current_setting.choice_hash[tag] = new_value
      return new_value
    end
    
    # choice of a boolean. default value is false
    def boolean(tag, default_value = false)
      choice(tag){|inst, val|
        if val.nil?
          default_value
        else
          STDERR.puts "used to be #{val}, norikae #{not val}"
          not val
        end
      }
    end
    
    # default value is top of the list
    def list(tag, candidate_list) 
      len = candidate_list.length
      throw "too short list for :#{tag}" if len <= 0
      choice(tag){|inst, val|
        if (val.nil?) || len == 1        
          candidate_list[0] 
        else
          ran = rand(len-1)
          if (ind = candidate_list.index(val)) && ran >= ind
            ran += 1
          end
          candidate_list[ran]
        end
      }
    end
    
    def positive_integer(tag, default_value)
      throw "default value nonpositive :#{default_value}" if default_value <= 0
      choice(tag){|inst, val|
        if val.nil?
          default_value
        else
          cand = val
          begin
            cand = (cand * Math::exp(2* initial_setting.evolution_instruction.drasticness*(2*rand()-1))).to_i
            cand = [1,cand].max
            cand += 2* rand(2)-1
            cand = [1,cand].max
          end while cand == val
          cand
        end
      }  
    end
    
    def positive_integer_multiple_of(tag, default_value, factor)
      return factor * positive_integer(tag, default_value / factor)
    end
    

    attr_accessor :initial_setting, :current_setting
  end
end
