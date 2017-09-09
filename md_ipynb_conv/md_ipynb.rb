# md_ipynb.rb
# convert md to ipynb with splitted by '#'

require 'pp'
require 'json'

source = File.read(ARGV[0])

cont = []
cont0 = ""
source.each_line do |line|
  if m=line.match(/^#/)
    cont << cont0
    cont << line
    cont0 = ""
  else
    cont0 += line
  end
end
cont << cont0

ipynb = JSON.load(File.read('template.ipynb'))

cont.each do |lines|
  md_cell0 = {"cell_type"=>"markdown", "metadata"=>{}, "source"=>[]}
  md_cell0["source"] << lines
  ipynb["cells"] << md_cell0
end

File.open('tmp.ipynb','w') do |target|
  JSON.dump(ipynb,target)
end
