def detect_machine ()
  fn = File.expand_path('~/.ssh/id_rsa.pub')
  
  unless File.exist?(fn)
    throw 'cannot find '+fn
  end
  
  key = open(fn,'r').read 

  return :nushiomac if key.index('nushio@muranushi-takashigyou-no-macbook.local')
  return :tsubame if key.index('10IKA101@t2a')
  return :degima if key.index('GDVSXkCaTlwVLsw')
  return :momigi if key.index('momigi')
  throw 'unknown machine'
end

Machine = detect_machine()
GPUPerNode = case Machine
             when :tsubame
               3
             when :momigi
               3
             when :degima
               2
             else
               1
             end
