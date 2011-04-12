#!/usr/bin/env ruby

require 'optparse'
require 'fileutils'
require 'nushio_script/optimization'

include Optimization

BathPath =  FileUtils.pwd + '/bath'
LogPath = BathPath + '/log'
ArchivePath = BathPath + '/archive'
BathListFn = BathPath + '/bath_list.bin'
ExecutableFn = 'auto_optimize.out'
BenchmarkDoneFn = 'benchmark.done'
BenchFn = 'bench.txt'


Option = {}

Nodes = (6161..6180).map{|i| "t2a00#{i}"}

$stderr = STDERR
$stdlog = open('/dev/null', 'w')
$free_nodes = []
$free_node_wait = 1
$free_node_limitter = 0

$most_popular_score = Unknown
$highscore = Unknown
$highscore_setting = nil

$bath_hash = {}
$bath_keys = []
def sh(str)
  $stderr.puts str
  system(str)
end

def get_free_node()
  if $free_nodes.length <= 0 && $free_node_limitter <= 0
    $free_node_limitter += 1
    $stderr.puts "polling for free nodes, sleep = #{$free_node_wait}..."
    sleep($free_node_wait)
    Nodes.each{|node|
      `ssh #{node} LANG=c uptime` =~ /load\s+average\:(.*)$/
      loads = $1.split(',').map{|x|x.to_f}
      if loads[0] < 1.0
        $free_nodes += [node]
      end
    }
    $free_nodes = $free_nodes.clone.sort_by{rand()}
    
    if $free_nodes.length <= 0
      $free_node_wait = [1, 2 * $free_node_wait].max
    else
      $free_node_wait = 1
    end
  end
  return nil if $free_nodes.length <= 0
  return  $free_nodes.shift
end

class Bath
  attr_accessor :key, :temperature, :evolution_instruction
  attr_accessor :last_setting, :current_setting
  attr_accessor :state, :idle_since, :running_since, :running_node
  attr_accessor :armageddon

  def initialize(key, temp, evolution_instruction) 
    @key = key
    @temperature = temp
    @evolution_instruction = evolution_instruction
    @last_setting = OptimizationSetting::new(@evolution_instruction)
    @current_setting = nil
    @state = :idle
    @idle_since = Time::now + rand()
    @running_since = Time::now + rand()
  end

  def tag()
    key.join('_')
  end
  
  def path()
    BathPath + '/' + tag()
  end
  
  def prepare()
    $stderr.puts "prepareing #{path()}..."
    `mkdir #{path()}`
    `svn export --force . #{path()}`
    `cp -r library/ #{path()}/library/` # have to materialize to avoid conflict
    `cp -r nushio_script/ #{path()}/nushio_script/`
    `ln -s $PWD/#{CanonFilename} #{path()}/#{CanonFilename}` if Option[:tsubame_gpu]
  end
  
  def take_track(msg)
    return unless Option[:keep_track]
    $stdlog.puts <<TRACK
#{Time::now().to_i} #{tag()} #{msg}
TRACK
  end

  def proceed()
    case @state
    when :idle
      if launch()
        @state = :benchmarking 
        take_track('launch')
        @running_since = Time::now()
      end
    when :benchmarking
      done_fn = path() + '/' + BenchmarkDoneFn
      unless File.exist?(done_fn) || (time_out = (Time::now() - @running_since > 300))
        @state = :benchmarking
      else
        take_track("timeout #{@running_node}") if time_out
        `rm -f #{done_fn}`
        `ssh #{@running_node} killall -9 #{ExecutableFn} &> /dev/null`
        open(path() + '/' + CurrentFilename) {|fp|
          @current_setting = Marshal.load(fp)          
        }
        score_num = 0.0; score_den = 0; error_count = 0;
        unless File.exist?(bench_fn=path() + '/' + BenchFn)
          $stderr.puts "no " + BenchFn
        else
          open(bench_fn, 'r') { |fp|
            while line = fp.gets
              judge, score = line.split(/\s+/)
              if judge != 'AC'
                error_count += 1
              else
                score_num += score.to_f
                score_den += 1
              end
            end
          }
        end
        final_score = 
          if score_den <= 0 || error_count > 0
            $stderr.puts "score_den <= 0" if score_den <= 0 
            $stderr.puts "error_count > 0" if error_count > 0
            0.0
          else
            score_num / score_den
          end
        @current_setting.score = final_score
        @state = :evaluated
        take_track("evaluated #{final_score}")
      end
    when :evaluated
      (@key.length).times{
        i = rand(@key.length)
        new_key = @key.clone
        new_key[i]+= 2*rand(2)-1
        unless (neighbour_bath = $bath_hash[new_key]).nil?
          if thermal_accept(neighbour_bath.last_setting)
            take_track("thermal_accept #{new_key.join('_')} #{neighbour_bath.last_setting.score}")
          end
        end
      }
      if thermal_accept(@current_setting)            
        take_track("thermal_accept self #{@current_setting.score}")
        if Option[:keep_archive] && @current_setting.score > $highscore
          arch_fn = "#{ArchivePath}/opt_#{tag()}_#{Time::now().to_i}.bin"
          open(arch_fn, 'w') {|fp| Marshal.dump(@current_setting, fp) }
        end
      end

      @state = :idle      
      @idle_since = Time::now
    else 
      $stderr.puts "bathtub is in unknown state :#{@state}"
    end
  end

  def thermal_accept(new_setting)
    return false if new_setting.score == Unknown
    if @last_setting.score==Unknown || 
        new_setting.score > @last_setting.score ||
        Math::exp( (new_setting.score -  @last_setting.score) / @temperature) > rand()
      @last_setting = new_setting
      return true
    end
    return false
  end

  def launch()
    node = get_free_node()
    
    return false unless node
    @running_node = node
    `rm -f #{path() + '/' + BenchFn}`

    batch_execute_fn = path() + '/batch_execute.sh'
    batch_launch_fn = path() + '/batch_launch.sh'
    
    fresh_setting = @last_setting.clone
    if @armageddon && fresh_setting.score == @armageddon
      take_track("armageddon #{$highscore}")
       @last_setting.score = fresh_setting.score = 0
    elsif fresh_setting.score == $most_popular_score && fresh_setting.score != Unknown
      $stderr.print "#{tag} have most popular score #{$most_popular_score} "
      if rand < 0.5
        $stderr.puts ": demo son na no kankei nee"
      elsif $most_popular_score == $highscore
        rand_key = $bath_keys[rand($bath_keys.length)]
        $stderr.puts ": merge with #{rand_key.join('_')}"
        fresh_setting.merge($bath_hash[rand_key].last_setting)
      else
        $stderr.puts ": merge with highscore"
        fresh_setting.merge($highscore_setting)
      end
    end
    open(path() + '/' + InitialFilename, 'w'){|fp|
      ts = fresh_setting.choice_hash.to_a.map{|t,v|t}
      henko_size = [1,(ts.length ** @evolution_instruction.broadness).to_i].max
      te = []
      henko_size.times{
        te << ts[rand(ts.length)]
      }
      fresh_setting.evolution_instruction = @evolution_instruction
      fresh_setting.tag_to_be_evolved = te
      Marshal.dump(fresh_setting,fp)
    }
    open(batch_execute_fn,'w'){ |fp|
      fp.puts <<SH
cd #{path()}
rake clobber
rake
./#{ExecutableFn} -g 0 -o #{BenchFn} &
sleep 1
./#{ExecutableFn} -g 1 -o #{BenchFn} &
sleep 1
./#{ExecutableFn} -g 2 -o #{BenchFn} 
touch #{BenchmarkDoneFn}
SH
    }
    log_instruction = if Option[:keep_log]
                        "&> #{LogPath}/#{tag()}_#{Time::now().to_i}.log"
                      else; '&> /dev/null'; end
    open(batch_launch_fn,'w'){ |fp|
      fp.puts <<SH
#{batch_execute_fn} #{log_instruction}
SH
    }
    `chmod 777 #{batch_execute_fn}`
    `chmod 777 #{batch_launch_fn}`
    $stderr.puts "submit job #{tag()} to node #{node}"
    system(<<SH)
ssh #{node} #{batch_launch_fn} &
SH
    return true
  end
  
end







parser = OptionParser.new
parser.on('-c', '--create', 'create optimization settings'){
  Option[:task] = :create
}
parser.on('-p', '--peek', 'prettyprint a opt_*.bin file'){
  Option[:task] = :peek
}
parser.on('-u', '--unmarshal', 'look inside of any Marshal-ed file.'){
  Option[:task] = :unmarshal
}
parser.on('-X', '--proceed', 'proceed 1 step of optimization'){
  Option[:task] = :proceed
}

parser.on('--track', 'keep track of happened events'){
  Option[:keep_track] = true
  $stdlog = open(BathPath + '/optimization.log', 'a')
}
parser.on('--log', 'keep every compile logs'){
  Option[:keep_log] = true
}
parser.on('--archive', 'keep every optimization parameters'){
  Option[:keep_archive] = true
}
parser.on('--stdout', 'output log to stdout instead of stderr'){
  $stderr = STDOUT
}
parser.on('--tsubame-gpu', 'set optimization environment for tsubame GPU') {
  Option[:tsubame_gpu] = true
}
parser.on('--unfreeze') {
  Option[:unfreeze] = true
}



parser.on('-t', '--temperature MAX MIN DIV', 'set bathtub temperatures. temperature list begins with MAX,' +
          ' each divided by DIV, while it is greater than MIN.') {|param|
  words = param.split(/\s+/)
  throw "specify three float for --temperature, e.g. --temperature '1.0 1e-3 3.16'" unless words.length == 3
  ma, mi, div = words.map{|x|x.to_f}
  t = ma; ret = []
  while t >= mi
    ret << t
    t /= div
  end
  
  Option[:temperature_list] = ret.sort
}


if ARGV.length <= 0
  $stderr.puts parser
  exit 0
end

parser.parse!(ARGV)

case Option[:task]
when :peek
  ARGV.each {|fn|
    open(fn, 'r') {|fp|
      set = Marshal.load(fp)
      params = set.choice_hash.to_a.sort.map{|tag, val| ":#{tag} = #{val}"}
      $stderr.puts <<SETTING
drasticness = #{set.evolution_instruction.drasticness}
broadness   = #{set.evolution_instruction.broadness}
#{params.join("\n")}
score       = #{set.score}
SETTING
    }
  }
when :unmarshal
  ARGV.each {|fn|
    open(fn, 'r') {|fp|
      puts Marshal.load(fp).inspect
    }
  }
when :create
  throw 'also need --temperature to create.' unless Option[:temperature_list]

  Nodes.each{|node|  sh "ssh #{node} killall -9 #{ExecutableFn} &"}
  `rm -fr #{BathPath}/*`
  `mkdir -p #{BathPath}`
  `mkdir -p #{LogPath}`
  `mkdir -p #{ArchivePath}`
  dra_list = Option[:drasticness_list] || [0, 0.5, 1.0] 
  bro_list = Option[:broadness_list] || [0, 0.25, 0.5]
  
  bro_list.length.times{|bro_id| bro = bro_list[bro_id]
    dra_list.length.times{|dra_id| dra = dra_list[dra_id]
      Option[:temperature_list].length.times{|temp_id| temp = Option[:temperature_list][temp_id]
        # the key is combination of three numbers
        # smaller first key means lower broadness.
        # smaller second key means lower drasticness.
        # smaller thired key means lower temperature.
 
        # in short, smaller keys indicate conservativeness / shortsight.
        
        key = [bro_id, dra_id, temp_id] 
        $bath_hash[key] = bath = Bath::new(key, temp, EvolutionInstruction::new(bro, dra))
        bath.prepare()
      }
    }
  }
  open(BathListFn,'w') {|fp|
    Marshal.dump($bath_hash.to_a.sort, fp)
  }
when :proceed
  $stderr.puts "Optimization in Progress " + Time::now.to_i.to_s
  score_histogram = {}
  open(BathListFn,'r') {|fp|
    bath_kvpair = Marshal.load(fp)
    bath_kvpair.each{|k,v| 
      if Option[:unfreeze]
        v.last_setting.score = Unknown
      end
      if Option[:temperature_list]
        v.temperature = Option[:temperature_list][v.key[2]]
      end
      }
    $bath_hash = {}
    bath_kvpair.each{|k,v| 
      $bath_hash[k] = v       
      $bath_keys.push(k)
      s = v.last_setting.score
      score_histogram[s] ||= 0
      score_histogram[s] += 1
    }
    $most_popular_score = score_histogram.to_a.map{|s, n| [n,-s]}.sort[-1][1]*(-1)
    $most_popularity = score_histogram.to_a.map{|s, n| [n,-s]}.sort[-1][0]
    $highscore = score_histogram.keys.sort[-1]
    bath_kvpair.each{|k,v| 
      $highscore_setting = v.last_setting.clone if  v.last_setting.score == $highscore
    }
  }


  $stderr.puts "most popular score is #{$most_popular_score}; popularity #{$most_popularity}/#{$bath_keys.length}; highscore is #{$highscore}"

  if 2*$most_popularity > $bath_keys.length && $most_popular_score == $highscore
    $stderr.puts "*** WRATH OF GOD ****"
    $bath_hash.each{|k,b|
      b.armageddon = $highscore
    }
  end

  $bath_hash.sort_by{|k,bath| bath.idle_since}.each{|key, bath|
    bath.proceed()
  }
  open(BathListFn,'w') {|fp|
    Marshal.dump($bath_hash.to_a.sort, fp)
  }
  
  if true
    msg_sum = ''
    $bath_hash.to_a.sort_by{|key,val| key.reverse}.each{|key, bath|
      msg = ''
      msg += bath.tag()
      #msg += " T=" + bath.temperature.to_s
      msg += " S=" + bath.last_setting.score.to_s[0..8] unless bath.last_setting.score == Unknown

      msg += " " + bath.state.to_s[0..2]
      msg += " " + bath.current_setting.score.to_s[0..8] if bath.state == :evaluated

      msg = msg[0...39]
      msg += ' ' while msg.length < 40
      msg_sum += msg
      msg_sum += "\n" if key[0]==2
    }
    $stderr.puts msg_sum
    open(BathPath+'/optimizer_state.txt','w') {|fp|  fp.puts msg_sum}
  end
  
end
