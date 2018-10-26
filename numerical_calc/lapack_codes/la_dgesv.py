import numpy as np
import scipy.linalg as la
import time

sizes = [1000, 2000, 4000]
for n in sizes:
    A = np.random.random((n, n))
    b = np.random.random((n))
    start = time.time()
    la.lapack.dgesv(A, b)
    print('%s [dim], %10.5f [sec] # python' % (n,time.time()-start))



