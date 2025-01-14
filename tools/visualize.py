import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread

original = mmread(sys.argv[1])
row_perm = np.loadtxt(sys.argv[2])
col_perm = np.loadtxt(sys.argv[3])
# row_perm = [ 6, 2, 4, 3, 5, 0, 1, 7 ]
# col_perm = [ 4, 1, 5, 2, 6, 3, 0, 7 ]

n = len(col_perm)
id_perm = range(0, n)

permutation_left = csr_matrix((np.ones(n, dtype=int), (row_perm, id_perm)), shape=(n, n))
permutation_right = csr_matrix((np.ones(n, dtype=int), (id_perm, col_perm)), shape=(n, n))

fig, axs = plt.subplots(1, 2)

axs[0].spy(original, markersize=0.1)
axs[1].spy(permutation_left @ original @ permutation_right, markersize=1)

plt.show()