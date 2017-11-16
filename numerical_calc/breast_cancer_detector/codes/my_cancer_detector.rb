require 'narray'
def read_n_data_from_lines(matrix_A, vector_b, lines)
  m,n =  matrix_A.shape
  lines[0..n-1].each_with_index do |line,i|
    row = line.chomp.split(',')
    vector_b[i] = (row[1] == 'B')? 1.0 : -1.0
    row[2..-1].each_with_index { |ele, j|
      matrix_A[j,i] = ele.to_f
    }
  end
end

def calc_dLw(matrix_A, vector_b, vector_w, opts)
  # calc L(w) from A,b and w
  m, n = matrix_A.shape
  vector_lw = NArray.sfloat(m)
  sum, correct = 0, 0
  n.times do |i|
    tmp = (matrix_A[0..-1,i] * vector_w)[0]
    tmp = (tmp<0)? -1.0 : 1.0 if opts[:step]
    correct += 1 if tmp*vector_b[i]>0
    tmp = tmp - vector_b[i]
    sum += tmp*tmp
  end

  print "L(w) = #{sum}\ncorrect = #{correct}/#{n}\n" if opts[:print]

  # calc dL(w)
  vector_Lw = NVector.sfloat(m)
  m.times do |j|
    n.times do |i|
      tmp = (matrix_A[0..-1,i] * vector_w)[0]
      tmp = (tmp<0)? -1.0 : 1.0 if opts[:step]
      tmp = tmp - vector_b[i]
      vector_Lw[j] += tmp*matrix_A[j,i]
    end
  end
  return vector_Lw
end 

def print_w(vector_w)
  params = ["radius", "texture", "perimeter","area","smoothness","compactness",
            "concavity","concave points","symmetry","fractal dimension"];
  print("      params      :      (mean)    (stderr)     (worst)\n")
  j = 0
  vector_w.each do |w_j|
    printf("%17s :",params[j/3]) if j%3==0
    printf("%12.4f", w_j)
    print "\n" if j%3==2
    j += 1
  end
end

# initial set ups
data_dir = '../Klein_codes'
lines = File.readlines(File.join(data_dir,'train.data'))

p n = lines.size
p m = lines[0].split(',').size - 2

matrix_A = NMatrix.sfloat(m,n)
vector_b = NVector.sfloat(n)
vector_w = NVector.sfloat(m)
vector_Lw = NVector.sfloat(m)

read_n_data_from_lines(matrix_A, vector_b, lines)

n.times do |i|
  m.times do |j|
    vector_Lw[j] += matrix_A[j,i]/n
  end
end
m.times do |j|
  vector_w[j] = 1.0
  #  vector_w[j] = 0.0
  # vector_w[j] = 1.0/vector_Lw[j]/30
end

print "# Finding minimum loop\n"
#loop, sigma = 300, 3.0*10**(-9)
loop, sigma = 100, 3.0*10**(-9)
printf("* Loop: %5d\n* sigma: %11.4e\n",loop,sigma)
step_sel = false
printf("* Use step function: %s\n\n", step_sel.to_s)
print "## Initial values of Vector w\n"
print_w(vector_w)
vector_Lw = calc_dLw(matrix_A, vector_b, vector_w, opts={:print => true, :step => step_sel})
print_w(vector_Lw)
loop.times do |l|
  print "." if l%10==0
  #  print "#{i}-loop\n"
  vector_w = vector_w - vector_Lw*sigma
  vector_Lw = calc_dLw(matrix_A, vector_b, vector_w, opts={:print => false, :step => step_sel})
end
print "\n"

print "## Adjusted values of Vector w\n"
print_w(vector_w)
vector_Lw = calc_dLw(matrix_A, vector_b, vector_w, opts={:print => true, :step => step_sel})

# validate data judgement
print "# Validate.data results\n"
lines = File.readlines(File.join(data_dir,'validate.data'))

n = lines.size
m = lines[0].split(',').size - 2

matrix_A = NMatrix.sfloat(m,n) #
vector_b = NVector.sfloat(n) # int

read_n_data_from_lines(matrix_A, vector_b, lines)

#print_w(vector_w)
vector_Lw = calc_dLw(matrix_A, vector_b, vector_w, opts={:print => true, :step => step_sel})
