module Config
  
  class Definition
    attr_accessor :key, :value
    def initialize(k,v)
      @key = k; @value = v
    end
    def to_s ()
      return :value
    end
  end
 
  def self.parse_config(str)
    def self.parse_macro(str)
      words = str.split(/\s+/)
      macro = words[0][1..-1]
      case macro
        when 'define' , 'undef'
        return [Definition::new(words[1], str)]
      else
        return [Definition::new('', str)]
      end
    end

    def self.parse_c(str0)
      def self.parse_single_c(str)
        return Definition::new('', str) if str.strip.length <= 0
        return Definition::new('', str) if str =~ /\}/
        return Definition::new('', str) if str =~ /\{/
        words = str.gsub(';',' ').gsub('=',' = ').split(/\s+/) 
        if words.index('typedef')
          return Definition::new(words[-1], str)
        elsif indeq = words.index('=')
          return Definition::new(words[indeq-1], str)
        else
          return Definition::new('', str)
        end
      end
      ret = []
      str0.split(/([^\;\}]*[\;\}])/).each{| str |
        ret << parse_single_c(str)
      }
      return ret
    end

    # recombinate macro
    recombinated = str.split(/\n/).map{|line|
      if line[-1..-1] == "\\"
        line
      else
        line + "\n"
      end
    }.join.split(/\n/).map{|x| x + "\n"}

    defs = []
    buf = ''
    recombinated.each{|line|
      if line[0..0] == '#'
        defs += parse_c(buf); buf = ''
        defs += parse_macro(line) 
      else
        buf += line 
      end
    }
    defs += parse_c(buf); buf = ''
    return defs
  end


  def self.generate(customization)
    defo_config = parse_config(open('library/simulation_config_default.h','r').read)
    cust_config = parse_config(customization)

    new_config = defo_config.clone
    
    cust_config.each{|my_def|
      next if my_def.key == ''
      ind = defo_config.index{ |your_def| 
        your_def.key == my_def.key
      }
      if not ind
        # STDERR.puts "Warning: new definition:  #{my_def.value}"
        new_config << my_def
      else
        new_config[ind].value = my_def.value
      end
    }
    
    return new_config.map{|x|
      x.value
    }.join
  end
  def self.generate!(x)
    puts generate(x)
  end
end
