#require 'narray'
require 'nmatrix/lapacke'
#require 'nmatrix'
# initial set ups
dir = '.' #'./codes'
lines_A = File.readlines(File.join(dir, 'train_A.data'))
lines_b = File.readlines(File.join(dir, 'train_b.data'))

p n = lines_A.size
p m = lines_A[0].split("\t").size
matrix_A = NMatrix.new([n,m], dtype: :float32,stype: :dense)
vector_b = NMatrix.new([n], dtype: :float32,stype: :dense)
vector_w = NMatrix.new([m], dtype: :float32,stype: :dense)

n.times do |i|
  lines_A[i].split("\t").each_with_index do |v,j|
    matrix_A[i,j] = v.to_f
  end
  vector_b[i] = lines_b[i].to_f
end

def print_w(vector_w)
  params = ["radius", "texture","perimeter","area",
            "smoothness","compactness","concavity","concave points",
            "symmetry","fractal dimension"];
  print("    (params)      :")
  print("    (mean)    (stderr)     (worst)")
  params.each_with_index do |param, i|
    printf("\n%17s :",param)
    3.times{|j| printf("%12.8f", vector_w[i*3+j])}
  end
end

qq, rr = matrix_A.factorize_qr
q = qq[0..n-1,0..m-1]
r = rr[0..m-1,0..m-1]

print(q.transpose.shape)
print(vector_b.shape)
q.transpose*vector_b

