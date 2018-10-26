require 'open3'


result = []

['1000','2000','4000'].each do |size|
  out, _err, _status = Open3.capture3("./lapack_argv #{size}")
  puts out
  m = out.match(/(\d+).+(\d+\.\d+)/)
  result << [m[1].to_i,m[2].to_f]
end

p result

