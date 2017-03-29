#!/opt/local/bin/python3.5

#'./dimReduce.py' NOT WORKING ON MY LOCAL COMPUTER

####################
#Reducing dimensions of full_matrix.txt
#
#Input: full_matrix.txt with full path
#
#Output: n x n principal component matrix
#that will be used for HMM learning
####################

from sklearn.decomposition import PCA
import numpy as np

#instantiate PCA object
file = '/Users/nickgiangreco/GitHub/HMMicro/data/final_matrix.txt'

matrix_full = np.loadtxt(file, dtype='i', delimiter='\t', skiprows = 1)
matrix_full_T = matrix_full.transpose
#print(matrix_full) #success!

pca = PCA(n_components=matrix_full.transpose().shape[1])
pca.fit(matrix_full.transpose())
matrix_reduce = pca.transform(matrix_full.transpose())

print(matrix_reduce.shape) #huzzah!

output = '/Users/nickgiangreco/GitHub/HMMicro/data/matrix_reduce.txt'
np.savetxt(output, matrix_reduce, delimiter="\t")



