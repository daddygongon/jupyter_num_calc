# files = ['introduction_scores_180619.txt','speaker_180619.org']
files = ['introduction_scores_180620.txt','speaker_180620.org']

report = File.readlines(files[0])
speaker = File.readlines(files[1])

speaker.each do |s_line|
  p s_line
  report.each_with_index do |r_line, i|
    if r_line.include?(s_line.chomp)
      p r_line
      r_split = r_line.chomp.split(',')
      r_split[0] = 10
      replace = r_split.join(',')
      report[i] = replace+"\n"
    end
  end
end

scores = []
report.each do |line|
  score = line.split(',')[0]
  line.split(',')[1..-1].each do |id|
    scores << [id,score]
  end
end

tmp = File.open("tmp.csv",'w')
luna_sheet = File.readlines('luna_sheet.tmp')
luna_sheet[1..-1].each do |line|
  warden = false
  scores.each do |score|
    if line.include?(score[0])
      tmp.print(line.chomp+","+score[0].chomp+","+score[1]+"\n")
      warden = true
    end
  end
  tmp.print(line) if warden==false
end
tmp.close
