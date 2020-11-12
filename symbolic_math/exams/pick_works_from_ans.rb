require 'json'

puts "ruby pick_workds_from_ans.rb FILE LIMIT=-1 MD_DELS=\'\' CODE_REMS=\'\'"
file = ARGV[0] || 'gw_1_exp_log_ans.ipynb'

limit = ARGV[1] || -1
limit = limit.to_i

dels = ARGV[2] || ''
p dels = dels.split(' ').map(&:to_i)

rems = ARGV[3] || ''
p rems = rems.split(' ').map(&:to_i)


p ['source: ', file]
json = JSON.parse(File.read(file))

# pick keys
json.each_pair do |key, cont|
  p key
end

# select delete
json["cells"][0..limit].each_with_index do |cont, i|
  p [i, cont["cell_type"], cont["source"][0]]
  next if rems.include?(i)
  dels << i if cont["cell_type"] == 'code'
end

p ['dels',dels.sort]
# delete code cell from json
dels.sort.reverse.each do |i|
  json["cells"].delete_at(i)
end

# dump file to 'gw_1_exp_log.ipynb'
dump_file = file.gsub('_ans.', '.')
p ['dump:  ', dump_file]
File.open(dump_file,'w'){|f| JSON.dump(json, f) }

