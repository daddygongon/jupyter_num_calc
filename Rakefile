puts 'help README.org and ipynbs, should copy this Rakefile to the cwd'

task :default do
  system 'rake -T'
end

desc 'sync FILE with hiki'
task :sync do
  p file = ARGV[1]
  entries = if file.nil?
              Dir.entries('.')[2..-1]
            else
              [file]
            end
  p attach_dir = '/Users/bob/Sites/new_ist_data/ist_data/cache/attach/DoingMathWithPython/'
  p http_address = 'https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math'
  entries.each do |file|
    ['.data','.pdf','.ipynb'].each do |ext|
      if File.extname(file)==ext
        print "{{attach_anchor(#{file})}}"
        print "([[nbviewer|#{http_address}/#{file}]])\n"
        FileUtils.cp(file, File.join(attach_dir,file))
      end
    end
  end
  exit
end

desc 'mk_link FILE for org'
task :mk_link do
  p file = ARGV[1]
  entries = if file.nil?
              Dir.entries('.')[2..-1]
            else
              [file]
            end
  p attach_dir = 'https://github.com/daddygongon/jupyter_num_calc/tree/master/symbolic_math'
  p http_address = 'https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math'
  entries.each do |file|
    #    ['.data','.ipynb'].each do |ext|
    ['.ipynb'].each do |ext|
      if File.extname(file)==ext
        print "- [[#{attach_dir}/#{file}][#{file}]]"
        print "([[#{http_address}/#{file}][nbviewer]])\n"
      end
    end
  end
  exit
end


desc 'conv_link FILE.org > FILE2.org'
task :conv_link do
  lines = File.readlines(ARGV[1])
# [[_d1][file]]([[address][_d2]])
lines.each do |line|
  if m=line.match(/\[\[(.*?)\]\[(.*?)\]\]\(\[\[(.*?)\]\[(.*?)\]\]\)/)
    _d1, file, address, _d2 = m[1..-1]
    puts  "- [[#{address}][#{file}]]"
  else
    puts line
  end
end

  dir = ARGV[1]
end

desc 'ls Dir'
task :ls do
  dir = ARGV[1]
  p Dir.entries(dir)
  exit
end
