"""
Performance evaluation measures for Protein Complex Prediction (PCP)
Author: Md. Shahidul Islam [shahidcseku@gmail.com]
Date: 12 July 2020
Description: Python implementation of performance evaluation measures for 
			Protein Complex Prediction (PCP) algorithms. 
			Measures included: Accuracy, Precision, Recall, 
			F-measure, Sensitivity, Positive Predictive Value (PPV), etc. 
			It can be also used for graph clustering.
"""

import math

def evaluatePerformance(output_file, benchmark_file, overlap, minClusterSize=2):
	# Parsing file
	output = open(output_file,"r")
	benchmark = open(benchmark_file,"r")
	output_data = output.read().splitlines()
	benchmark_data = benchmark.read().splitlines()

	# build clusters
	clusters = 0
	output_clusters = []
	total_predicted_proteins = 0
	for output_cluster in output_data:
		string_list = output_cluster.split('\t')
		string = [x for x in string_list if not x.isdigit()]
		total_predicted_proteins += len(string)

		if(len(string)>minClusterSize):
			output_clusters.append(string)
			clusters += 1
		# endif
	# endfor

	benchmark_clusters = []
	benchmark_proteins = 0

	for benchmark_cluster in benchmark_data:
		string_list = benchmark_cluster.split('\t')
		string = [x for x in string_list if not x.isdigit()]
		benchmark_proteins += len(string)
		# Not necessary, in case needed
		if(len(string)>minClusterSize):
			benchmark_clusters.append(string)
		# endif
		benchmark_clusters.append(string)
	# endfor

	# Find metric components
	P_abs = len(output_clusters)
	B_abs = len(benchmark_clusters)
	Ncp = 0
	Ncb = 0
	for index, output_cluster in enumerate(output_clusters):
		for index2, benchmark_cluster in enumerate(benchmark_clusters):
			affinity = NeighborhoodAffinity(output_cluster, benchmark_cluster)
			if(affinity>=overlap):
				Ncp += 1
				break
			# endif
		# endfor
	# endfor

	for index, benchmark_cluster in enumerate(benchmark_clusters):
		for output_cluster in (output_clusters):
			affinity = NeighborhoodAffinity(output_cluster, benchmark_cluster)
			if(affinity>=overlap):
				Ncb += 1
				break
			# endif
		# endfor
	# endfor

	# For SEN
	maxiTij = [0 for x in range(len(benchmark_clusters))]
	Bi_sum = 0
	for index, benchmark_cluster in enumerate(benchmark_clusters):
		Bi_sum += len(benchmark_cluster)
		for index2, output_cluster in enumerate(output_clusters):
			score = matchingScore(output_cluster, benchmark_cluster)
			if(maxiTij[index]<score):
				maxiTij[index] = score
			# endif
		# endfor
	# endfor

	# For PPV
	maxjTij = [0 for x in range(len(output_clusters))]
	Oi_sum = 0
	U_BiPj = 0
	for index, output_cluster in enumerate(output_clusters):
		Oi_sum += len(output_cluster)
		U = []
		for index2, benchmark_cluster in enumerate(benchmark_clusters):
			score = matchingScore(output_cluster, benchmark_cluster)
			if(maxjTij[index]<score):
				maxjTij[index] = score
			# endif
			U = set(U).union(set(output_cluster) & set(benchmark_cluster))
		# endfor
		U_BiPj += len(U)
	# endfor

	# print measures

	print("P", P_abs)
	print("B", B_abs)
	print('avg cluster size', total_predicted_proteins/P_abs)
	print("# Clusters", clusters)
	print("Ncp", Ncp)
	print("Ncb", Ncb)
	print("Covered proteins", sum(maxiTij))

	Sensitivity = sum(maxiTij)/Bi_sum
	PPV = sum(maxjTij)/(U_BiPj)
	Accuracy = math.sqrt(Sensitivity*PPV)

	print("Sensitivity ", Sensitivity)
	print("PPV ", PPV)
	print("Accuracy ", Accuracy)

	Precesion = Ncp/P_abs
	Recall = Ncb/B_abs
	F_measure = (2*Precesion*Recall)/(Precesion+Recall)

	print("Precesion ", Precesion)
	print("Recall ", Recall)
	print("F-measure ", F_measure)
# endfunction

# similarity funciton
def NeighborhoodAffinity(p, b):
	score = matchingScore(p, b)
	affinity = (score*score)/(len(p)*len(b))
	return affinity
# endfunction

def matchingScore(cluster, benchmark, overlap=0.2):
	score = len(set(benchmark) & set(cluster))
	return score
# endfunction


# RUN =========================================================
# Adjust your settings
overlap = 0.20
minClusterSize = 3

# Change filename with directory
filename = "input/sample.txt"
benchmark = "benchmark/cmplx.3.cyc2008.tab.txt"


print("Performance Evaluation for Protein Complex Prediction")
print("=====================================================")
evaluatePerformance(filename,benchmark, overlap, minClusterSize)