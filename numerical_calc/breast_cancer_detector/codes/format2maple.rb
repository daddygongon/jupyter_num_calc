data_dir = '.'
source = 'validate' #'train'
lines = File.readlines(File.join(data_dir,"#{source}.data"))

target = File.open("#{source}_A.data",'w')
lines.each do |line|
  data = line.split(',')
  data[2..-2].each {|data| target.printf("%f\t",data)}
  target.printf("%f\n",data[-1])
end
target.close

target = File.open("#{source}_b.data",'w')
lines.each do |line|
  data = line.split(',')
  res = (data[1]=='M')? -1 : 1
  target.puts res
end
target.close
