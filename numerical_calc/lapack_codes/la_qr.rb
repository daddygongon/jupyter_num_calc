require 'nmatrix/lapacke'
#require 'nmatrix'
require 'time'

[1000,2000,4000].each do |size|
  n = size
  aa = NMatrix.random([n,n+1])
  start = Time.now
  #  lu = aa.factorize_lu
  q,r  = aa.factorize_qr
 
#  NMatrix::LAPACK.gesvd(aa)
  p Time.now - start
end
