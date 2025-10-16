file = 'wdbc.data'
data = File.readlines(file).map { |line| line.chomp.split(',') }


[[0, 299, 'train'], [300, data.size - 1, 'validate']
].each do |init, finish, label|
  b_data = []
  a_data = []
  data[init..finish].each do |row|
    if row[1] == 'B'
      b_data << -1 # 良性 benign
    else
      b_data << 1 # 悪性 malignant
    end
    a_data << row[2..-1].map(&:to_f).join(' ')
  end

  File.write("#{label}_b.data", b_data.join("\n"))
  File.write("#{label}_A.data", a_data.join("\n"))
end
