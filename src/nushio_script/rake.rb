require 'rake/clean'

task :default => :build
desc 'Build programs'
task :build 

desc 'Perform tests on updated parts only'
task :test 
desc 'Perform tests on entire tests'
task :test_all 
desc 'Perform tests, ignoring all errors'
task :test_force



CXX = 'nvcc'
NVCC = 'nvcc'
LINKER = 'nvcc'

SRC_EXTS = ['cpp', 'c', 'cu', 'inl', 'h']
OBJ_EXT = 'o'
EXE_EXT = 'out'

# cprb option.
# note that this default rule assumes 
# that cprb generates the same content for all instances of the library.
# if you need problem-specific generation of Source files,
# you may manually specify so.
SRC_EXTS.each{|exten|

  `find -L . -name '*.#{exten}rb'`.split(/\n/).each{|fn0| 
    fn = fn0[2..-1]
    target = fn.ext(exten)
    file target => fn do |t|
      sh "rm -rf #{t.name}"
      sh "cprb #{fn0} --keep -o #{t.name}"
      sh "chmod 444 #{t.name}"
    end
    CLEAN.include(target)
  }
}

class RakeOption
  attr_accessor :compile, :link, :compiler, :linker, :obj_tag
  def initialize()
    @compile = ['-I.', '-I$HOME/.nvcc/include']   # compiler options
    @link = []      # linker options
    @compiler = nil # force compiler
    @linker = nil   # force linker
    @obj_tag = []   # object file tags
  end
end

def recurse(dir, tasks)
  unless File.exist?("#{dir}/Rakefile")
    STDERR.puts "Warning: no Rakefile: #{dir}"
    return
  end
  tasks.each{|task_tag|
    task task_tag do
      sh("cd #{dir}; rake -N #{task_tag}")
    end
  }
end

# specify inclusion notation
def cinclude(hash)
  file hash do
    # cause the file to be 'updated' automatically
    # when the included files are updated
    sh "touch #{hash.to_a[0][0]}"
  end
end

# register tasks to compile one, or multiple files.
# returns the list of the resultant object files
def compile(fns, option = RakeOption::new)
  ret = []
  [fns].flatten.each {|fn|
    obj_fn = (option.obj_tag+[fn.ext('')]).join('-').gsub('/','-').ext(OBJ_EXT)
    case File.extname(fn)
    when '.cpp', '.cxx', '.cc'
      file obj_fn => fn do |t|
        sh "#{option.compiler || CXX} #{fn} -c -o #{obj_fn} #{option.compile.join(' ')}"
      end
    when '.cu'
      file obj_fn => fn do |t|
        sh "#{option.compiler || NVCC} #{fn} --keep -c -o #{obj_fn} #{option.compile.join(' ')}"
        # leave .cubin and .ptx , remove any other intermediate files
        im_files = ['cudafe?.*', 'linkinfo', 'cu.cpp', 'fatbin.*', 'hash', 'cpp?.*'].map{|excess_ext|
          fn.ext(excess_ext)
        }
        sh "rm -f #{im_files.join(' ')}"
      end
      CLEAN.include(fn.ext('ptx'))
      CLEAN.include(fn.ext('*.cubin'))
    end
    CLEAN.include(obj_fn)
    ret.push(obj_fn)
  }
  return ret
end

# register tasks to link req and create target
def link(target, req, option = RakeOption::new)
  file target => req do |task|
    sh "#{option.linker || CXX} #{task.prerequisites.join(' ')} -o #{task.name} #{option.link.join(' ')}"
  end
  CLOBBER.include(target)
  return target
end

GTEST_COMPILE_FLAGS = '-I'+`gtest-config --includedir`.strip
GTEST_LINK_FLAGS = '-L'+`gtest-config --libdir`.strip+' -lgtest -lgtest_main ' +
                   `gtest-config --ldflags`.strip
# register test suite
def gtest(test_name, req_cpps, option0 = RakeOption::new)
  target_fn = test_name.ext(EXE_EXT)
  option = option0.clone
  option.obj_tag.push(test_name)
  option.compile.push(GTEST_COMPILE_FLAGS)
  option.link.push(GTEST_LINK_FLAGS)
  req_objs = compile(req_cpps, option)
  link(target_fn, req_objs, option)
  
  passport_fn = File.dirname(target_fn) + '/.passed.' + File.basename(target_fn)

  task :build => target_fn
  task :test => [passport_fn, target_fn]
  file passport_fn => target_fn do
    sh "./#{target_fn}"
    sh "touch #{passport_fn}"
    CLEAN.include(passport_fn)
  end
  task :test_all => target_fn do
     sh "./#{target_fn}"
  end
  task :test_force => target_fn do
    ret = sh "./#{target_fn}" do  |ok, res|
    end
  end
end

# register vanilla c-like program
def cprogram(program_name, req_cpp_objs, option0 = RakeOption::new)
  option = option0.clone
  option.obj_tag.push(program_name)
  target_fn = program_name.ext(EXE_EXT)
  req_objs = []
  req_cpps = []
  req_cpp_objs.each{|fn|
    if SRC_EXTS.index(File.extname(fn)[1..-1])
      req_cpps.push(fn)
    elsif OBJ_EXT == File.extname(fn)[1..-1]
      req_objs.push(fn)
    else
      STDERR.puts "Unknown ext: #{fn}"
      exit(-1)
    end
  }
  req_objs += compile(req_cpps, option)
  link(target_fn, req_objs, option)
  task :build => target_fn
end
