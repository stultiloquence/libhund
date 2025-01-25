import json
import pprint
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread

def load_bisection_result(matrix, parallel_attempts, max_imbalance):
	filename = f"../testing/results/bisection-test_{matrix}_{parallel_attempts}_{str(max_imbalance).rjust(2, "0")}.json"
	print(filename + '\n')
	with open(filename) as f: result = json.load(f)
	return result

matrices = [ "raefsky6", "lns_3937", "fd18", "twotone" ]
# matrices = [ "fd18" ]
parallel_attempts_all = [ 64 ]
max_imbalances = [ 5, 10, 20, 40 ]
results = {}
for matrix in matrices:
	results[matrix] = {}
	for parallel_attempts in parallel_attempts_all:
		results[matrix][parallel_attempts] = {}
		for max_imbalance in max_imbalances:
			results[matrix][parallel_attempts][max_imbalance] = load_bisection_result(matrix, parallel_attempts, max_imbalance)


for matrix in matrices:
	y = [ min(results[matrix][64][imbalance][0]["relative_separator_sizes"]) for imbalance in max_imbalances ]
	print(matrix, y)
	plt.plot(max_imbalances, y, label=matrix)
plt.xscale("log")
plt.legend()
plt.show()