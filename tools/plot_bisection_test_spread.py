import json
import pprint
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmread

def load_bisection_result(matrix, parallel_attempts, max_imbalance):
	filename = f"../testing/spread/bisection-test_{matrix}_{parallel_attempts}_{str(max_imbalance).rjust(2, "0")}.json"
	print(filename + '\n')
	with open(filename) as f: result = json.load(f)
	return result

matrices = [ "raefsky6", "lns_3937", "fd18", "twotone", "cage14" ]
parallel_attempts_all = [ 32 ]
max_imbalances = [ 5 ]
results = {}
for matrix in matrices:
	results[matrix] = {}
	for parallel_attempts in parallel_attempts_all:
		results[matrix][parallel_attempts] = {}
		for max_imbalance in max_imbalances:
			results[matrix][parallel_attempts][max_imbalance] = load_bisection_result(matrix, parallel_attempts, max_imbalance)


_, axes = plt.subplots(1, len(matrices))

matplotlib.rcParams.update({'font.size': 22})
matplotlib.rc('ytick', labelsize=20)

for i, matrix in enumerate(matrices):
	# for imbalance in max_imbalances:
		# print(results[matrix][64][imbalance][0]["relative_separator_sizes"])
	values = results[matrix][32][5][0]["relative_separator_sizes"]
	values.sort()
	axes[i].scatter([0] * len(values), values, label=matrix, s=100)
	axes[i].get_xaxis().set_visible(False)
	axes[i].legend()
	axes[i].yaxis.set_tick_params(labelsize=20)
	# minimum = values[0]
	# maximum = values[-1]
	# median = values[len(values) // 2]
	# ratio = minimum / median
	# print(f"Matrix {matrix} has min {minimum}, max {maximum}, median {median}, min/median = {ratio}.")

plt.legend()
plt.show()