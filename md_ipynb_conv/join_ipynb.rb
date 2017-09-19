# join_ipynb.rb
# join ipynbs 

require 'pp'
require 'json'

ipynb0 = JSON.load(File.read(ARGV[0]))
ipynb1 = JSON.load(File.read(ARGV[1]))

ipynb0["cells"].each do |cell|
  pp cell
  ipynb1["cells"] << cell
end

File.open('tmp.ipynb','w') do |target|
  JSON.dump(ipynb1,target)
end
