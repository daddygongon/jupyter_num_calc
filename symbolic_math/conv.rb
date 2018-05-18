lines = File.readlines(ARGV[0])

# [[https://github.com/daddygongon/jupyter_num_calc/tree/master/symbolic_math/cg.ipynb][cg.ipynb]]([[https://nbviewer.jupyter.org/github/daddygongon/jupyter_num_calc/blob/master/symbolic_math/cg.ipynb][nbviewer]])
lines.each do |line|
  if m=line.match(/\[\[(.*?)\]\[(.*?)\]\]\(\[\[(.*?)\]\[(.*?)\]\]\)/)
    _d1, file, address, _d2 = m[1..-1]
    puts  "- [[#{address}][#{file}]]"
  else
    puts line
  end
end
